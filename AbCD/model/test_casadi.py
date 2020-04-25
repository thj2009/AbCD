#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
from casadi import *
import numpy as NP
import matplotlib.pyplot as plt
from operator import itemgetter

nk = 20  # Control discretization
tf = 10.0  # End time

# Declare variables (use scalar graph)
u = SX.sym("u")  # control
x = SX.sym("x", 2)  # states

# ODE right hand side and quadratures
xdot = vertcat([(1 - x[1] * x[1]) * x[0] - x[1] + u, x[0]])
qdot = x[0] * x[0] + x[1] * x[1] + u * u

# DAE residual function
dae = SXFunction("dae", daeIn(x=x, p=u), daeOut(ode=xdot, quad=qdot))

# Create an integrator
integrator = Integrator("integrator", "cvodes", dae, {"tf": tf / nk})

# All controls (use matrix graph)
x = MX.sym("x", nk)  # nk-by-1 symbolic variable
U = vertsplit(x)  # cheaper than x[0], x[1], ...

# The initial state (x_0=0, x_1=1)
X = MX([0, 1])

# Objective function
f = 0

# Build a graph of integrator calls
for k in range(nk):
    X, QF = itemgetter('xf', 'qf')(integrator({'x0': X, 'p': U[k]}))
    f += QF

# Terminal constraints: x_0(T)=x_1(T)=0
g = X

# Allocate an NLP solver
nlp = MXFunction("nlp", nlpIn(x=x), nlpOut(f=f, g=g))
solver = NlpSolver("solver", "ipopt", nlp)

# Solve the problem
sol = solver({"lbx": -0.75,
              "ubx": 1,
              "x0": 0,
              "lbg": 0,
              "ubg": 0})

# Retrieve the solution
u_opt = NP.array(sol["x"])

# Time grid
tgrid_x = NP.linspace(0, 10, nk + 1)
tgrid_u = NP.linspace(0, 10, nk)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid_u, u_opt, '-.')
plt.title("Van der Pol optimization - single shooting")
plt.xlabel('time')
plt.legend(['u trajectory'])
plt.grid()
plt.show()

import sys
import tempfile

import numpy as np
import casadi as cas
from .network import KineticModel
from AbCD.utils import get_index_species


class TPDcondition(object):
    def __init__(self, name=''):
        self.name = name
        self.T0 = 100
        self.Tf = None
        self.Beta = 10
        self.SimulationTime = None
        self.Ntime = 1000
        self.TimeGrid = None
    def _calSim(self):
        self.SimulationTime = (self.Tf - self.T0) / self.Beta
    def _calTimeGrid(self):
        self.TimeGrid = np.linspace(0, self.SimulationTime, self.Ntime)

class VacTPDcondition(TPDcondition):
    def __init__(self, name=''):
        TPDcondition.__init__(self, name)
        self.InitCoverage = {}

class VacuumTPD(KineticModel):
    '''
    Vacuum Temperature Programmed Desorption
    '''
    def __init__(self, specieslist, reactionlist, dEa_index, dBE_index):
        KineticModel.__init__(self, specieslist, reactionlist, dEa_index, dBE_index)
        self.pressure_value = []
        self.coverage_value = []
        self.rate_value = {}
        self.equil_rate_const_value = {}
        self.energy_value = {}

        self._t = cas.SX.sym('t', 1)                # Time
        self._T0 = cas.SX.sym('T0', 1)              # Initial Temperature
        self._beta = cas.SX.sym('beta', 1)          # Beta: Temperature Changing Rate
        # Derivative and Algebraic variable
        self._d_partP = cas.SX.sym('d_partP', self.ngas)
        self._d_cover = cas.SX.sym('d_cover', self.nsurf)

    def initialize(self, scale=1.0, pump_level=1):
        '''
        Build the dae system of transient CSTR model
        :param tau: space time of CSTR model
        '''
        self._Tem = self._T0 + self._beta * self._t
        self.build_kinetic(constTem=500)
        self.build_rate(scale=1)
        for j in range(self.ngas):
            self._d_partP[j] = - pump_level * self._partP[j]
            for i in range(self.nrxn):
                self._d_partP[j] += - 1 / pump_level * self.stoimat[i][j] * self._rate[i]
        for j in range(self.nsurf):
            self._d_cover[j] = 0
            for i in range(self.nrxn):
                self._d_cover[j] += self.stoimat[i][j + self.ngas] * self._rate[i]
        self._x = cas.vertcat([self._partP, self._cover])
        self._p = cas.vertcat([self._dEa, self._dBE, self._T0, self._beta])
        self._xdot = cas.vertcat([self._d_partP, self._d_cover])
#        self._dae_ = dict(x=self._x, p=self._p, ode=self._xdot, t=self._t)
        self._dae_ = cas.SXFunction("dae", cas.daeIn(x=self._x,p=self._p,t=self._t),
                                    cas.daeOut(ode=self._xdot))

    def fwd_simulation(self, dE_start, condition, detail=True,
                       reltol=1e-6, abstol=1e-8):

        time = condition.TimeGrid
        T0 = condition.T0
        beta = condition.Beta

        opts = {}
        opts['abstol'] = abstol
        opts['reltol'] = reltol
        opts['disable_internal_warnings'] = True
        opts['max_num_steps'] = 1e5

        x0 = self.init_condition(condition)
        P_dae = np.hstack([dE_start, T0, beta])
#        print(x0)
#        print(P_dae)
#        print(time)
        opts['tf'] = 2
        Fint = cas.Integrator('Fint', 'cvodes', self._dae_, opts)
        F_sim = Fint(x0=x0, p=P_dae)


#        Fint = cas.Integrator('Fint', 'cvodes', self._dae_, opts)
#        Fsim = cas.Simulator('Fsim', Fint, time)
#        Fsim.setInput(x0, 'x0')
#        Fsim.setInput(P_dae, 'p')
#        Fsim.evaluate()
        
#        pressure = Fsim.getOutput().full()

    def init_condition(self, condition):
        x0 = [0] * (self.nspe - 1)
        for spe in condition.InitCoverage.keys():
            idx = get_index_species(spe, self.specieslist)
            x0[idx] = condition.InitCoverage[spe]
        x0.append(1 - sum(x0))
        return x0
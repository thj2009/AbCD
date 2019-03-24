import sys
import tempfile

import numpy as np
import casadi as cas
from .network import KineticModel
from AbCD.utils import get_index_species


class BATCHcondition(object):
    def __init__(self, name=''):
        self.name = name
        self.Temperature = 298.15       # unit Kelvin
        self.TotalPressure = 1          # total pressure, unit: atm
        self.TotalSite = 0.00           # unit: total mole of activesite
        self.SimulationTime = 0
        self.Ntime = 2000
        self.PartialPressure = {}        # Dictionary to store the partial pressure of gas phase species
        self.TurnOverFrequency = {}      # Dictionary to store outlet turnover frequency
        self.InitRate = {}
    def _calGrid(self):
        self.TimeGrid = np.logspace(-5, self.SimulationTime, self.Ntime, base=10)

class Batch(KineticModel):
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
        # Derivative and Algebraic variable
        self._d_partP = cas.SX.sym('d_partP', self.ngas)
        self._d_cover = cas.SX.sym('d_cover', self.nsurf)

    def initialize(self, scale=1.0, A_V_ratio=1e-6, constTem=None):
        '''
        Build the dae system of transient CSTR model
        :param tau: space time of CSTR model
        '''
        self.build_kinetic(constTem=constTem)
        self.build_rate(scale=1)
        for j in range(self.ngas):
            self._d_partP[j] = 0
            for i in range(self.nrxn):
                self._d_partP[j] += A_V_ratio * self.stoimat[i][j] * self._rate[i]

        for j in range(self.nsurf):
            self._d_cover[j] = 0
            for i in range(self.nrxn):
                self._d_cover[j] += self.stoimat[i][j + self.ngas] * self._rate[i]
        self._x = cas.vertcat([self._partP, self._cover])
        self._p = cas.vertcat([self._dEa, self._dBE, self._Tem])
        self._xdot = cas.vertcat([self._d_partP, self._d_cover])
        self._dae_ = cas.SXFunction("dae", cas.daeIn(x=self._x, p=self._p, t=self._t),
                                    cas.daeOut(ode=self._xdot))

    def fwd_simulation(self, dE_start, condition, detail=True,
                       reltol=1e-6, abstol=1e-8):

        Tem = condition.Temperature
        time = condition.TimeGrid

        opts = {}
        opts['abstol'] = abstol
        opts['reltol'] = reltol
        opts['disable_internal_warnings'] = True
        opts['max_num_steps'] = 1e5

        P_dae = np.hstack([dE_start, Tem])


        # Partial Pressure
        Pinlet = np.zeros(self.ngas)
        for idx, spe in enumerate(self.specieslist):
            if spe.phase == 'gaseous':
                Pinlet[idx] = condition.PartialPressure[str(spe)]
        x0 = Pinlet.tolist() + [0] * (self.nsurf - 1) + [1]
        # print(x0)


        Fint = cas.Integrator('Fint', 'cvodes', self._dae_, opts)
        Fsim = cas.Simulator('Fsim', Fint, time)
        Fsim.setInput(x0, 'x0')
        Fsim.setInput(P_dae, 'p')
        Fsim.evaluate()
        
        # Evalu
        out = Fsim.getOutput().full()
        return out
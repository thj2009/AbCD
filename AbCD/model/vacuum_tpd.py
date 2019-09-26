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
        self.TemGrid = None
        self.TemProfile = None
        self.RateProfile = None
        self.PeakPosition = None

    def __repr__(self):
        return self.name

    def _calSim(self):
        self.SimulationTime = (self.Tf - self.T0) / self.Beta
    def _calGrid(self):
        self.TimeGrid = np.linspace(0, self.SimulationTime, self.Ntime)
        self.TemGrid = np.linspace(self.T0, self.Tf, self.Ntime)

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

        self._prior_ = None
        self.prob_func = None

    def initialize(self, scale=1.0, pump_level=1e5, A_V_ratio=1e-3, des_scale=1e-3, constTem=None):
        '''
        Build the dae system of transient CSTR model
        :param tau: space time of CSTR model
        '''
        self._Tem = self._T0 + self._beta * self._t
        self.build_kinetic(constTem=constTem)
        self.build_rate(scale=1, des_scale=des_scale)
        for j in range(self.ngas):
            self._d_partP[j] = - pump_level * self._partP[j]
            for i in range(self.nrxn):
                self._d_partP[j] += A_V_ratio * self.stoimat[i][j] * self._rate[i]
        for j in range(self.nsurf):
            self._d_cover[j] = 0
            for i in range(self.nrxn):
                self._d_cover[j] += self.stoimat[i][j + self.ngas] * self._rate[i]
        self._x = cas.vertcat([self._partP, self._cover])
        self._p = cas.vertcat([self._dEa, self._dBE, self._T0, self._beta])
        self._xdot = cas.vertcat([self._d_partP, self._d_cover])
        self._dae_ = cas.SXFunction("dae", cas.daeIn(x=self._x, p=self._p, t=self._t),
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
#        opts['tf'] = 2
#        Fint = cas.Integrator('Fint', 'cvodes', self._dae_, opts)
#        F_sim = Fint(x0=x0, p=P_dae)


        Fint = cas.Integrator('Fint', 'cvodes', self._dae_, opts)
        Fsim = cas.Simulator('Fsim', Fint, time)
        Fsim.setInput(x0, 'x0')
        Fsim.setInput(P_dae, 'p')
        Fsim.evaluate()
        
        # Evaluate
        out = Fsim.getOutput().full()
        return out

    def init_condition(self, condition):
        x0 = [0] * (self.nspe - 1)
        for spe in condition.InitCoverage.keys():
            idx = get_index_species(spe, self.specieslist)
            x0[idx] = condition.InitCoverage[spe]
        x0.append(1 - sum(x0[self.ngas:]))
        return x0


    def eval_prob(self, dE, conditionlist, evidence_info, prior_info):
        # evaluate prior
        self.prior_construct(prior_info)
        prob_func = self.prob_func
        prob_func.setInput(dE, 'i0')
        prob_func.evaluate()
        log_prior = -float(prob_func.getOutput('o0'))
        # evaluate likihood
        log_likeli = self.eval_likeli(dE, conditionlist, evidence_info)
        return log_likeli, log_prior

    def eval_likeli(self, dE, conditionlist, evidence_info={}):
        reltol = evidence_info.get('reltol', 1e-12)
        abstol = evidence_info.get('abstol', 1e-12)

        err = evidence_info.get('peak_err', 10)

        opts = {}
        opts['abstol'] = abstol
        opts['reltol'] = reltol
        opts['disable_internal_warnings'] = True
        opts['max_num_steps'] = 1e5

        # Initialize simulator
        evidence = 0
        for condition in conditionlist:
            time = condition.TimeGrid
            T0 = condition.T0
            beta = condition.Beta
            x0 = self.init_condition(condition)

            P_dae = np.hstack([dE, T0, beta])

            Fint = cas.Integrator('Fint', 'cvodes', self._dae_, opts)
            Fsim = cas.Simulator('Fsim', Fint, time)
            Fsim.setInput(x0, 'x0')
            Fsim.setInput(P_dae, 'p')
            Fsim.evaluate()

            out = Fsim.getOutput().full()
            # Find the peak
            for spe, peak_exp in condition.PeakPosition.items():
                idx = get_index_species(spe, self.specieslist)
                des = out[idx, :]
                idx_peak = np.argmax(des)
                peak_sim = condition.TemGrid[idx_peak]
                dev = peak_sim - peak_exp
                evidence += (dev * dev)/err**2
        return -evidence

    def prior_construct(self, prior_info):
        Pnlp = self._Pnlp
        if prior_info['type'] == 'Ridge':
            L2 = prior_info['L2']
            prior = cas.mul(Pnlp.T, Pnlp) * L2
        elif prior_info['type'] == 'Gaussian':
            mean = prior_info['mean']
            cov = prior_info['cov']
            dev = Pnlp - mean
            prior = cas.mul(cas.mul(dev.T, np.linalg.inv(cov)), dev)
        elif prior_info['type'] == 'GP':
            BEmean = prior_info['BEmean']
            BEcov = prior_info['BEcov']
            linear_BE2Ea = prior_info['BE2Ea']
            Eacov = prior_info['Eacov']

            _BE = Pnlp[self._NEa:]
            _Ea = Pnlp[:self._NEa]
            Eamean = cas.mul(linear_BE2Ea, _BE)
            dev_BE = _BE - BEmean
            dev_Ea = _Ea - Eamean
            prior = cas.mul(cas.mul(dev_BE.T, np.linalg.inv(BEcov)), dev_BE) + \
                    cas.mul(cas.mul(dev_Ea.T, np.linalg.inv(Eacov)), dev_Ea)
        self._prior_ = prior
        self.prob_func = cas.MXFunction('prob_func', [Pnlp], [prior])
        return prior
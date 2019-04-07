import sys
import tempfile

import numpy as np
import casadi as cas
from .network import KineticModel

class DynamicCSTRCondition(object):
    '''
    CSTR reaction condition class
    '''
    def __init__(self, name=''):
        self.name = name
        self.Temperature = 298.15
        self.TotalFlow = 0.00           # unit: s-1
        self.TotalSite = 0.00           # unit: total mole of activesite
        self.SimulationTime = 0
        self.Ntime = 1000
        self.NtimeSS = 100
        self.timeSS = 10000
        self.PartialPressure = {}        # Dictionary to store the partial pressure of gas phase species

    def __repr__(self):
        return self.name

    def detail_repr(self):
        '''
        Return the detial representation of the reaction condition
        '''
        out = ''
        out += self.name
        out += 'Temperature = %.2f K' %(self.Temperature)
        out += 'Total Flow rate = %.2e s-1' %(self.TotalFlow)
        out += 'Simulation time = %.2e s' %(self.SimulationTime)
        out += 'Inlet Pressure (atm): %.2e' %(self.TotalPressure)
        for key, value in self.PartialPressure.items():
            out += '    {0:<10s} {1:<10.2e}'.format(key, value)
        out += 'Experiment Outlet'
        for key, value in self.TurnOverFrequency.items():
            out += '    {0:<10s} {1:<10.2e}'.format(key, value)
        out += '\n'
        return out


class DynamicCSTR(KineticModel):
    '''
    CSTR modeling TOOL
    '''
    def __init__(self, specieslist, reactionlist, dEa_index, dBE_index):
        KineticModel.__init__(self, specieslist, reactionlist, dEa_index, dBE_index)
        self.pressure_value = []
        self.coverage_value = []
        self.rate_value = {}
        self.equil_rate_const_value = {}
        self.energy_value = {}

        self.timeGrid = None
        self.output = None
        self.control = None

        self._t = cas.SX.sym('t', 1)                # Time
        self._partP_in = cas.SX.sym('partP_in', self.ngas)
        self._Flowtot = cas.SX.sym('Flowtot', 1)
        self._flow = cas.SX.sym('flow', self.ngas)
        # Derivative and Algebraic variable
        self._d_flow = cas.SX.sym('d_flow', self.ngas)
        self._d_cover = cas.SX.sym('d_cover', self.nsurf)

    def initialize(self, tau=1.0, scale=1.0):
        '''
        Build the dae system of transient CSTR model
        :param tau: space time of CSTR model
        '''
        Ptot = cas.sumRows(self._partP_in)            # Total Pressure
        for i in range(self.ngas):
            self._partP[i] = self._flow[i]/self._Flowtot * Ptot
        self.build_kinetic()
        self.build_rate(scale=1)
        for j in range(self.ngas):
            self._d_flow[j] = self._partP_in[j] / Ptot * self._Flowtot \
                            - self._flow[j]
            for i in range(self.nrxn):
                self._d_flow[j] += self.stoimat[i][j] * self._rate[i]
            # Scale the flow rate
            self._d_flow[j] /= tau
        for j in range(self.nsurf):
            self._d_cover[j] = 0
            for i in range(self.nrxn):
                self._d_cover[j] += self.stoimat[i][j + self.ngas] * self._rate[i]
        self._x = cas.vertcat([self._flow, self._cover])
        # ???
        self._p = cas.vertcat([self._dEa, self._dBE])
        self._u = cas.vertcat([self._partP_in, self._Tem, self._Flowtot])
        self._xdot = cas.vertcat([self._d_flow, self._d_cover])
        #self._dae_ = dict(x=self._x, p=self._p, ode=self._xdot)
        self._dae_ = cas.SXFunction('dae', cas.controldaeIn(x=self._x, p=self._p, u=self._u, t=self._t),
                                   cas.daeOut(ode=self._xdot))

    def fwd_simulation(self, dE_start, condition, detail=True,
                       reltol=1e-8, abstol=1e-10):

        # Simulation output grid
        N = condition.Ntime
        timeSS = condition.timeSS
        NSS = condition.NtimeSS
        if condition.timeSS != 0:
            tgrid = np.append(np.linspace(0, timeSS-1, NSS),
                              np.linspace(timeSS, timeSS + condition.SimulationTime, N))
        else:
            tgrid = np.linspace(0, condition.SimulationTime, N)
#        print(tgrid.shape)
        # ControlSimulator will output on each node of the timegrid
        opts = {}
        opts["integrator"] = "cvodes"
        opts["integrator_options"] = {"abstol": abstol,
                                      "reltol": reltol,
                                      'disable_internal_warnings': True,
                                      "max_num_steps": 1e5}
        Fsim = cas.ControlSimulator("Fsim", self._dae_, tgrid, opts)
#        print(Fsim.getInput("u").shape)

        # Initialize Coverage
        x0 = [0] * (self.nspe - 1) + [1]
        # Initialize Parameter
        P_dae = np.hstack([dE_start])
        # Initialize Control Variable
        control = self._calControledVar(condition)

        Fsim.setInput(x0, "x0")
        Fsim.setInput(P_dae, "p")
        Fsim.setInput(control, "u")
        Fsim.evaluate()

        self.timeGrid = tgrid[NSS:]
        self.control = control[:, NSS:]
        self.output = Fsim.getOutput().full()[:, NSS:]

    def _calControledVar(self, condition):
        '''
        y = A sin(B(x + C)) + D
        '''
        N = condition.Ntime - 1
        NSS = condition.NtimeSS
        tgrid = np.linspace(0, condition.SimulationTime, N)
        
        ssPinlet =np.zeros([self.ngas, NSS])
        Pinlet = np.zeros([self.ngas, N])
        for idx, spe in enumerate(self.specieslist):
            if spe.phase == 'gaseous':
                oscillateType = condition.PartialPressure[str(spe)]['type']
                partPressure = condition.PartialPressure[str(spe)]['partialPressure']
                ssPinlet[idx, :] = partPressure
                if oscillateType == 'constant':
                    Pinlet[idx, :] = partPressure
                elif oscillateType == 'sine':
                    period = condition.PartialPressure[str(spe)]['period']
                    phase = condition.PartialPressure[str(spe)]['phase']
                    amplitude = condition.PartialPressure[str(spe)]['amplitude']
                    Pinlet[idx, :] = partPressure + amplitude * np.sin((2 * np.pi) / period * (tgrid + phase)) 
                else:
                    pass
        # Dynamic State
        Tem = np.ones(N) * condition.Temperature
        totalFlow = np.ones(N) * condition.TotalFlow
        ctrl = np.vstack([Pinlet, Tem, totalFlow])
        # Steady State
        ssTem = np.ones(NSS) * condition.Temperature
        ssTotalFlow = np.ones(NSS) * condition.TotalFlow 
        ssCtrl = np.vstack([ssPinlet, ssTem, ssTotalFlow])
        ctrl = np.hstack([ssCtrl, ctrl])
        return ctrl
        
        
        

        # result = {}
        # # Detailed Reaction network data
        # # Evaluate partial pressure and surface coverage
        # self.pressure_value = (F_sim['xf'][:self.ngas]/TotalFlow * TotalPressure).full().T[0].tolist()
        # self.coverage_value = F_sim['xf'][self.ngas:].full().T[0].tolist()
        
        
        # # Evaluate Reaction Rate automatically save to Rate attribute
        # x = self._x
        # p = self._p
        # rate_fxn = cas.SXFunction('rate_fxn', [x, p],
                                  # [self._rate,
                                  # self._rfor,
                                  # self._rrev])
        # rate_fxn.setInput(F_sim['xf'], 'i0')
        # rate_fxn.setInput(P_dae, 'i1')
        # rate_fxn.evaluate()
        # self.rate_value['rnet'] = rate_fxn.getOutput('o0').full().T[0].tolist()
        # self.rate_value['rfor'] = rate_fxn.getOutput('o1').full().T[0].tolist()
        # self.rate_value['rrev'] = rate_fxn.getOutput('o2').full().T[0].tolist()
        
        # # Evaluate Reaction Energy
        # ene_fxn = cas.SXFunction('ene_fxn', [x, p],
                                # [self._reaction_energy_expression['activation'],
                                # self._reaction_energy_expression['enthalpy']])
        # ene_fxn.setInput(F_sim['xf'], 'i0')
        # ene_fxn.setInput(P_dae, 'i1')
        # ene_fxn.evaluate()
        # self.energy_value['activation'] = ene_fxn.getOutput('o0').full().T[0].tolist()
        # self.energy_value['enthalpy'] = ene_fxn.getOutput('o1').full().T[0].tolist()
        
        # # Evaluate Equilibrium Constant and Rate Constant
        # k_fxn = cas.SXFunction('k_fxn', [x, p],
                               # [self._Keq,
                               # self._Qeq,
                               # self._kf,
                               # self._kr])
        # k_fxn.setInput(F_sim['xf'], 'i0')
        # k_fxn.setInput(P_dae, 'i1')
        # k_fxn.evaluate()
        # self.equil_rate_const_value['Keq'] = k_fxn.getOutput('o0').full().T[0].tolist()
        # self.equil_rate_const_value['Qeq'] = k_fxn.getOutput('o1').full().T[0].tolist()
        # self.equil_rate_const_value['kf'] = k_fxn.getOutput('o2').full().T[0].tolist()
        # self.equil_rate_const_value['kr'] = k_fxn.getOutput('o3').full().T[0].tolist()
        
        # # RESULT
        # result['pressure'] = self.pressure_value
        # result['coverage'] = self.coverage_value
        # result['rate'] = self.rate_value
        # result['energy'] = self.energy_value
        # result['equil_rate'] = self.equil_rate_const_value

        #return tor, result
    
    # def condilist_fwd_simul(self, dE_start, conditionlist,
                            # reltol=1e-10, abstol=1e-12):
        # #Simulate the condition list and return result
        # tor_list = []
        # result_list = []
        
        # for condi in conditionlist:
            # tor, result = self.fwd_simulation(dE_start, condi, reltol=1e-10, abstol=1e-12)
            # tor_list.append(tor)
            # result_list.append(result)
        # return tor_list, result_list



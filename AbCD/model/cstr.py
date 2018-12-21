import sys
import tempfile

import numpy as np
import casadi as cas
from .network import KineticModel

class CSTRCondition(object):
    '''
    CSTR reaction condition class
    '''
    def __init__(self, name=''):
        self.name = name
        self.Temperature = 298.15       # unit Kelvin
        self.TotalPressure = 1          # total pressure, unit: atm
        self.TotalFlow = 0.00           # unit: s-1
        self.TotalSite = 0.00           # unit: total mole of activesite
        self.SimulationTime = 0
        self.PartialPressure = {}        # Dictionary to store the partial pressure of gas phase species
        self.TurnOverFrequency = {}      # Dictionary to store outlet turnover frequency

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


class CSTR(KineticModel):
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

    def initialize(self, scale=1.0, tau=1.0):
        self.build_kinetic()
        self.build_rate(scale=scale)
        self.build_dae(tau=tau)

    def build_dae(self, tau=1.0):
        '''
        Build the dae system of transient CSTR model
        :param tau: space time of CSTR model
        '''
        Ptot = cas.sumRows(self._partP_in)            # Total Pressure
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
        self._p = cas.vertcat([self._dEa, self._dBE, self._partP_in, self._Tem, self._Flowtot])
        self._xdot = cas.vertcat([self._d_flow, self._d_cover])
        self._dae_ = dict(x=self._x, p=self._p, ode=self._xdot)

    def fwd_simulation(self, dE_start, condition, detail=True,
                       reltol=1e-8, abstol=1e-10):

        TotalPressure = condition.TotalPressure
        TotalFlow = condition.TotalFlow
        Tem = condition.Temperature
        tf = condition.SimulationTime

        opts = {}
        opts['tf'] = tf       # Simulation time
        opts['abstol'] = abstol
        opts['reltol'] = reltol
        opts['disable_internal_warnings'] = True
        opts['max_num_steps'] = 1e5

        Fint = cas.Integrator('Fint', 'cvodes', self._dae_, opts)

        x0 = [0] * (self.nspe - 1) + [1]

        # Partial Pressure
        Pinlet = np.zeros(self.ngas)
        for idx, spe in enumerate(self.specieslist):
            if spe.phase == 'gaseous':
                Pinlet[idx] = condition.PartialPressure[str(spe)]
        P_dae = np.hstack([dE_start, Pinlet, Tem, TotalFlow])
        F_sim = Fint(x0=x0, p=P_dae)

        tor = {}
        for idx, spe in enumerate(self.specieslist):
            if spe.phase == 'gaseous':
                tor[str(spe)] = float(Pinlet[idx]/TotalPressure * TotalFlow - F_sim['xf'][idx])

        result = {}
        # Detailed Reaction network data
        # Evaluate partial pressure and surface coverage
        self.pressure_value = (F_sim['xf'][:self.ngas]/TotalFlow * TotalPressure).full().T[0].tolist()
        self.coverage_value = F_sim['xf'][self.ngas:].full().T[0].tolist()
        
        
        # Evaluate Reaction Rate automatically save to Rate attribute
        x = self._x
        p = self._p
        rate_fxn = cas.SXFunction('rate_fxn', [x, p],
                                  [self._rate,
                                  self._rfor,
                                  self._rrev])
        rate_fxn.setInput(F_sim['xf'], 'i0')
        rate_fxn.setInput(P_dae, 'i1')
        rate_fxn.evaluate()
        self.rate_value['rnet'] = rate_fxn.getOutput('o0').full().T[0].tolist()
        self.rate_value['rfor'] = rate_fxn.getOutput('o1').full().T[0].tolist()
        self.rate_value['rrev'] = rate_fxn.getOutput('o2').full().T[0].tolist()
        
        # Evaluate Reaction Energy
        ene_fxn = cas.SXFunction('ene_fxn', [x, p],
                                [self._reaction_energy_expression['activation'],
                                self._reaction_energy_expression['enthalpy']])
        ene_fxn.setInput(F_sim['xf'], 'i0')
        ene_fxn.setInput(P_dae, 'i1')
        ene_fxn.evaluate()
        self.energy_value['activation'] = ene_fxn.getOutput('o0').full().T[0].tolist()
        self.energy_value['enthalpy'] = ene_fxn.getOutput('o1').full().T[0].tolist()
        
        # Evaluate Equilibrium Constant and Rate Constant
        k_fxn = cas.SXFunction('k_fxn', [x, p],
                               [self._Keq,
                               self._Qeq,
                               self._kf,
                               self._kr])
        k_fxn.setInput(F_sim['xf'], 'i0')
        k_fxn.setInput(P_dae, 'i1')
        k_fxn.evaluate()
        self.equil_rate_const_value['Keq'] = k_fxn.getOutput('o0').full().T[0].tolist()
        self.equil_rate_const_value['Qeq'] = k_fxn.getOutput('o1').full().T[0].tolist()
        self.equil_rate_const_value['kf'] = k_fxn.getOutput('o2').full().T[0].tolist()
        self.equil_rate_const_value['kr'] = k_fxn.getOutput('o3').full().T[0].tolist()
        
        # RESULT
        result['pressure'] = self.pressure_value
        result['coverage'] = self.coverage_value
        result['rate'] = self.rate_value
        result['energy'] = self.energy_value
        result['equil_rate'] = self.equil_rate_const_value

        return tor, result
    
    def condilist_fwd_simul(self, dE_start, conditionlist,
                            reltol=1e-10, abstol=1e-12):
        # Simulate the condition list and return result
        tor_list = []
        result_list = []
        
        for condi in conditionlist:
            tor, result = self.fwd_simulation(dE_start, condi, reltol=1e-10, abstol=1e-12)
            tor_list.append(tor)
            result_list.append(result)
        return tor_list, result_list
            
    def evidence_construct(self, dE_start, conditionlist, evidence_info, res_rsample=False):

        err_type = evidence_info['type']
        err = evidence_info['err']

        Pnlp = self._Pnlp
        opts = fwd_sensitivity_option()
        Fint = cas.Integrator('Fint', 'cvodes', self._dae_, opts)

        x0 = [0] * (self.nspe - 1) + [1]
        evidence = 0
        for condition in conditionlist:
            TotalPressure = condition.TotalPressure
            TotalFlow = condition.TotalFlow
            Tem = condition.Temperature

            # Partial Pressure
            Pinlet = np.zeros(self.ngas)
            for idx, spe in enumerate(self.specieslist):
                if spe.phase == 'gaseous':
                    Pinlet[idx] = condition.PartialPressure[str(spe)]
            P_dae = cas.vertcat([Pnlp, Pinlet, Tem, TotalFlow])
            F_sim = Fint(x0=x0, p=P_dae)
            for idx, spe in enumerate(self.specieslist):
                if spe.phase == 'gaseous':
                    # TOR
                    tor = Pinlet[idx]/TotalPressure * TotalFlow - F_sim['xf'][idx]
                    if res_rsample:
                        rsample = np.random.normal(1,0.1,1)
                    else:
                        rsample = 1
                    if str(spe) in condition.TurnOverFrequency.keys():
                        if err_type == 'abs':
                            dev = tor - condition.TurnOverFrequency[str(spe)] * rsample
                        elif err_type == 'rel':
                            dev = 1 - tor/(condition.TurnOverFrequency[str(spe)] * rsample)
                        evidence += (dev * dev)/err**2
        self._evidence_ = evidence
        return evidence

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
        return prior

    def mle_estimation(self, dE_start, conditionlist, evidence_info, prior_info,
                       res_rsample=False,
                       nlptol=1e-2, maxiter=500, bfgs=True, print_level=5,
                       print_screen=False, report=''):
        Pnlp = self._Pnlp
        # Objective
        obj = self.evidence_construct(dE_start, conditionlist, evidence_info, res_rsample) +\
            self.prior_construct(prior_info)

        nlp = dict(f=obj, x=Pnlp, g=self._thermo_constraint_expression)

        nlpopts = {}
        nlpopts['max_iter'] = maxiter
        nlpopts['tol'] = nlptol
        nlpopts['acceptable_tol'] = 1e-3
        nlpopts['jac_d_constant'] = 'yes'
        nlpopts['expect_infeasible_problem'] = 'yes'
        nlpopts['hessian_approximation'] = 'limited-memory'
        nlpopts['print_level'] = print_level

        solver = cas.NlpSolver('solver', 'ipopt', nlp, nlpopts)

        # FIXIT
        lbP = np.array(prior_info['lbound'])
        ubP = np.array(prior_info['ubound'])
        # Thermo dynamic consistency check
        lbG = 0 * np.ones(2 * self.nrxn)
        ubG = np.inf * np.ones(2 * self.nrxn)

        x0 = np.hstack([dE_start])
        solution = solver(x0=x0, lbx=lbP, ubx=ubP, lbg=lbG, ubg=ubG)
        opt_sol = solution['x'].full().T[0].tolist()
        obj = solution['f'].full()[0][0]

        print('=' * 20)
        print('Starting Point:')
        print(dE_start)
        print('Parameter:')
        print(opt_sol)
        print('Objective:')
        print(obj)
        return opt_sol, obj

    def bayesian_infer(self, ntot, nbuf, dE_start, transi_matrix, 
                       conditionlist, evidence_info, prior_info,
                       sample_method='elementwise',
                       step_write=100, save_result=True):
        
        out = ''
        Pnlp = self._Pnlp
        # Objective
        likeli = self.evidence_construct(dE_start, conditionlist, evidence_info)
        prior = self.prior_construct(prior_info)
        
        prob_fxn = cas.MXFunction('prob_fxn', [Pnlp], [likeli, prior])
        prob_fxn.setInput(dE_start, 'i0')
        prob_fxn.evaluate()
        
        likeli_prev = np.exp(-float(prob_fxn.getOutput('o0')))
        prior_prev = np.exp(-float(prob_fxn.getOutput('o1')))
        posterior_prev = likeli_prev * prior_prev
        
        print('{0:^10s}  {1:^15s}  {2:^15s}  {3:^15s}  {4:^15s}'.\
              format('step', 'prior', 'likelihood', 'posterior', 'accept%'))
        print('==' * 30)
        
        tor_dis, result_dis = [], []
        if save_result:
            tor_prev, result_prev = self.condilist_fwd_simul(dE_start, conditionlist)

        Edis = []
        Eprev = np.copy(dE_start)
        
        jump = 0
        for i in range(ntot):
            if i % step_write == 0:
                print('{0:^10d}  {1:^15.2e}  {2:^15.2e}  {3:^15.2e}  {4:^15.2f}'.\
                        format(i, prior_prev, likeli_prev, posterior_prev,
                               jump/float(i+1)*100))
            Estar = _sample(Eprev, transi_matrix, sample_method)
            
            if not self.CheckThermoConsis(Estar) and\
                len(np.argwhere(np.array(Estar)-np.array(prior_info['lbound']) < 0)) != 0 and\
                len(np.argwhere(np.array(Estar)-np.array(prior_info['ubound']) > 0)) != 0:
                # REJECT
                pass
            else:
                try:
                    # Evaluate the posterior distribution
                    prob_fxn.setInput(Estar, 'i0')
                    prob_fxn.evaluate()
                    
                    likeli_star = np.exp(-float(prob_fxn.getOutput('o0')))
                    prior_star = np.exp(-float(prob_fxn.getOutput('o1')))
                    posterior_star = likeli_star * prior_star
                    if save_result:
                        for condi in conditionlist:
                            pass
                    # Determine accept or reject
                    UU = np.random.uniform()
                    AA = min(1, posterior_star/posterior_prev)
                    if UU < AA:
                        # ACCEPT
                        Eprev = np.copy(Estar)
                        prior_prev = prior_star
                        likeli_prev = likeli_star
                        posterior_prev = posterior_star
                        # OPTION: Evaluate the reaction kinetics
                        if save_result:
                            tor_prev, result_prev = self.condilist_fwd_simul(Estar, conditionlist)
                        jump+=1
                    else:
                        # REJECT
                        pass
                except:
                    print('not feasible')
                    pass
            if i >= nbuf:
                Edis.append(Eprev)
                if save_result:
                    tor_dis.append(tor_prev)
                    result_dis.append(result_prev)
        
        print('==' * 30)
        return Edis, tor_dis, result_dis
        
        
def _sample(dE_start, transi_matrix, sample_method):
    if sample_method == 'elementwise':
        n = len(dE_start)
        i = int(n * np.random.uniform())
        scale = np.random.normal(0,1)
        deltaE = scale * transi_matrix[:,i]
    if sample_method == 'augment':
        n = len(dE_start)
        scale = np.random.normal(0, 1, n)
        deltaE = transi_matrix.dot(scale)
    newE = np.array(dE_start) + np.array(deltaE)
    return list(newE)
            
        
        
        
def fwd_sensitivity_option(tf=5000, abstol=1e-10, reltol=1e-8):
    '''
    Options pass to CVODES for sensitivity analysis
    '''
    opts = {}
    opts['tf'] = tf
#    opts["linear_solver"] = "csparse"
    #opts["linear_solver_type"] = "user_defined"
    #opts['t0'] = 0
    #opts['print_stats'] = True
    opts['fsens_all_at_once'] = True
    opts['fsens_err_con'] = True
    opts['fsens_abstol'] = 1e-6
    opts['fsens_reltol'] = 1e-5
    opts['abstol'] = abstol
    opts['reltol'] = reltol
    opts['abstolB'] = 1e-6
    opts['reltolB'] = 1e-4
    #opts['ad_weight'] = 0
#    opts['ad_weight_sp'] = 1
#    opts['linear_multistep_method'] = 'bdf'
#    opts['exact_jacobian'] = False
#    opts['linear_solverB'] = 'lapacklu'
    #opts['linear_solver_typeB'] = 'dense'
#    opts['iterative_solverB'] = 'gmres'
    #opts['interpolation_type'] ='polynomial'
    #opts['max_multistep_order'] = 8
    opts['use_preconditioner'] = True
    opts['use_preconditionerB'] = True
    opts['pretype'] = 'both'
    opts['pretypeB'] = 'both'
    opts['steps_per_checkpoint'] = 1e3
    #opts['nonlinear_solver_iteration'] = 'functional'
    #opts['linear_solver_typeB'] = 'iterative' 
    opts['disable_internal_warnings'] = True
    #opts['sensitivity_method'] = 'staggered'
#    opts['max_num_steps'] = 1e5
    opts['stop_at_end'] = True
    return opts

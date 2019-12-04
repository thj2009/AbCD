# Reaction Network Class
import sys
import tempfile
import numpy as np
import casadi as cas

from AbCD.utils import get_index_species, check_index, find_index
from AbCD.utils import Constant as _const

class ReactionNet(object):
    '''
    Basic reaction network class
    '''
    def __init__(self, specieslist, reactionlist):
        self.specieslist = specieslist
        self.reactionlist = reactionlist
        self.nspe = len(specieslist)
        self.nrxn = len(reactionlist)
        self.ngas = len([spe for spe in self.specieslist if spe.phase == 'gaseous'])
        self.nsurf = len([spe for spe in self.specieslist if spe.phase == 'surface'])
        self.stoimat = self._stoimat_generate()

    def __repr__(self):
        rxnnet = ''
        for rxn in self.reactionlist:
            rxnnet += str(rxn) + '\n'
        return rxnnet

    def _stoimat_generate(self):
        stoimat = np.zeros([self.nrxn, self.nspe], dtype='int8')
        for i in range(self.nrxn):
            for react in self.reactionlist[i].reactant:
                k = get_index_species(str(react[0]), self.specieslist)
                stoimat[i, k] = -react[1]
            for prod in self.reactionlist[i].product:
                k = get_index_species(str(prod[0]), self.specieslist)
                stoimat[i, k] = prod[1]
        return stoimat

    def detail_repr(self):
        '''
        Detailed Reaction Network Representation
        Summary on Reactant and Product, and matrix like stoichiometry
        '''
        out = ''
        out += 'Reaction Network Summary' + '\n'
        out += 15 * self.nspe * '-' + '\n'
        out += 'Number of Species   : %d' %self.nspe + '\n'
        out += 'Number of Reactions : %d' %self.nrxn + '\n'
        out += 15 * self.nspe * '-' + '\n'
        out += '{:25s}'.format(' ')
        for spe in self.specieslist:
            out += '{:^10s}'.format(spe)
        out += '\n'
        for j, rxn in enumerate(self.reactionlist):
            out += '{:<25s}'.format(rxn)
            for i, spe in enumerate(self.specieslist):
                out += '{:^10d}'.format(self.stoimat[j, i])
            out += '\n'
        return out

    def plot_network(self):
        '''
        Plot the reaction network using pygraphviz
        '''
        # TODO:
        pass

class SimpleKinetic(ReactionNet):
    '''
    Simple Kinetic Class where the species energy and reaction kinetic is calcualted
    and stored.
    '''
    def __init__(self, specieslist, reactionlist, dEa_index, dBE_index):
        ReactionNet.__init__(self, specieslist, reactionlist)
        self.dEa_index = dEa_index
        self.dBE_index = dBE_index
        self._NEa = len(self.dEa_index)
        self._NBE = len(self.dBE_index)
        self._Np = self._NEa + self._NBE

        # Parameter
        self._dEa = cas.SX.sym('dEa', self._NEa)
        self._dBE = cas.SX.sym('dBE', self._NBE)
        self._Pnlp = cas.MX.sym('Pnlp', self._Np)
        # Controlled Variable
        self._Tem = cas.SX.sym('Tem', 1)
        # State variable
        self._cover = cas.SX.sym('cover', self.nsurf)
        self._partP = cas.SX.sym('partP', self.ngas)    # Partial Pressure of gas species
        # Kinetic Intermediate term, may used for functon evaluation
        self._thermo_constraint_expression = None
        self._reaction_energy_expression = None
        self._kf = None
        self._kr = None
        self._Keq = None
        self._Qeq = None
        self._dH_expression = None
        self._rate = None
        self._rfor = None
        self._rrev = None

    def build_kinetic(self, constTem=None, thermoTem=298.15):
        if constTem is not None:
            Tem = constTem
        else:
            Tem = self._Tem
        enthal_spe, enthal_spe_cons = [], []
        entro_spe = []
        denthal = []    # Store the enthalpy of each species for function calculation
        for j, spe in enumerate(self.specieslist):
            HH = spe.Enthalpy(Tem)
            SS = spe.Entropy(Tem)
            enthal_spe.append(HH)
            enthal_spe_cons.append(spe.Enthalpy(thermoTem))
            entro_spe.append(SS)
            denthal.append(0)
        # surface species add deviation variable
        for j, dbe_idx in enumerate(self.dBE_index):
            enthal_spe[dbe_idx] += self._dBE[j]
            enthal_spe_cons[dbe_idx] += self._Pnlp[self._NEa + j]
            denthal[dbe_idx] += self._Pnlp[self._NEa + j]

        kf = cas.SX.sym('kf', self.nrxn)
        kr = cas.SX.sym('kr', self.nrxn)
        Keq = cas.SX.sym('Keq', self.nrxn)
        enthal_react, entro_react, Ea_react = [], [], []
        enthal_react_cons, Ea_cons = [], []          # Store the Enthalpy of each reaction
        denthal_react = []
        for i, rxn in enumerate(self.reactionlist):
            kine_data = rxn.Arrhenius(Tem)
            Hreact, dHreact = 0, 0
            Sreact = 0
            Hreact_cons = 0
            Arr = kine_data['A']
            Ea0 = kine_data['Ea']
            Ea298 = rxn.Arrhenius(thermoTem)['Ea']
            n = kine_data['n']
            if check_index(i, self.dEa_index):
                ind = find_index(i, self.dEa_index)
                Ea = Ea0 + self._dEa[ind]
                Ea_cons.append(Ea298 + self._Pnlp[ind])
            else:
                Ea = Ea0
                Ea_cons.append(Ea298)
            Ea_react.append(Ea)

            kf[i] = Arr * (self._Tem**n) * \
                cas.exp(-Ea * 1000 / (_const.Rg * self._Tem))  # k = A * T^n * exp(-Ea/(RT))
            for j in range(self.nspe):
                Hreact += self.stoimat[i][j] * enthal_spe[j]
                dHreact += self.stoimat[i][j] * denthal[j]
                Hreact_cons += self.stoimat[i][j] * enthal_spe_cons[j]
                Sreact += self.stoimat[i][j] * entro_spe[j]

            enthal_react.append(Hreact)
            denthal_react.append(dHreact)
            enthal_react_cons.append(Hreact_cons)
            entro_react.append(Sreact)

            Keq[i] = cas.exp(-Hreact * 1000 / (_const.Rg * self._Tem)) * \
                    cas.exp(Sreact / _const.Rg)
            kr[i] = kf[i]/Keq[i]

        # Constraint function on parameter
        ineq2 = cas.vertcat(Ea_cons) - cas.vertcat(enthal_react_cons)
        thermal_consis_ineq = cas.vertcat([cas.vertcat(Ea_cons), ineq2])
        self._thermo_constraint_expression = thermal_consis_ineq
        self._reaction_energy_expression = \
            {'enthalpy': cas.vertcat(enthal_react), 'activation': cas.vertcat(Ea_react)}
        self._kf = kf
        self._kr = kr
        self._Keq = Keq
        self._dH_expression = denthal_react

    def build_rate(self, scale=1.0, des_scale=1):
        '''
        build net, forward, reverse rate expression for DAE system
        :param: scale: scale the rate constant for degree of rate control calculation
        '''
        # Net, Forward, Reverse rate for each elementary step
        rate = cas.SX.sym('rate', self.nrxn)         # s-1
        rfor = cas.SX.sym('rfor', self.nrxn)         # s-1
        rrev = cas.SX.sym('rrev', self.nrxn)         # s-1
        Qeq = cas.SX.sym('Qeq', self.nrxn)

        for i in range(self.nrxn):
            # Scale the rate constant
            rfor[i] = self._kf[i] * scale
            rrev[i] = self._kr[i] * scale
            Qeq[i] = 1
            for j in range(self.nspe):
                if self.stoimat[i][j] < 0:
                    if j < self.ngas:
                        partP = self._partP[j] * des_scale
                        rfor[i] *= partP**(-self.stoimat[i][j])
                        Qeq[i] /= partP**(-self.stoimat[i][j])
                    else:
                        rfor[i] *= (self._cover[j - self.ngas])**(-self.stoimat[i][j])
                        Qeq[i] /= (self._cover[j - self.ngas])**(-self.stoimat[i][j])
                elif self.stoimat[i][j] > 0:
                    if j < self.ngas:
                        partP = self._partP[j] * des_scale
                        rrev[i] *= partP**(self.stoimat[i][j])
                        Qeq[i] *= partP**(self.stoimat[i][j])
                    else:
                        rrev[i] *= (self._cover[j - self.ngas])**(self.stoimat[i][j])
                        Qeq[i] *= (self._cover[j - self.ngas])**(self.stoimat[i][j])
            rate[i] = rfor[i] - rrev[i]
        self._Qeq = Qeq
        self._rate = rate
        self._rfor = rfor
        self._rrev = rrev

    def CheckThermoConsis(self, dE_start, tol=1e-8):
        '''
        True: satify the thermodynamic consistency
        False: not satisfy
        '''
        Pnlp = self._Pnlp
        thermo_consis_fxn = cas.MXFunction('ThermoConsisFxn', [Pnlp], [self._thermo_constraint_expression])
        thermo_consis_fxn.setInput(dE_start, 'i0')
        thermo_consis_fxn.evaluate()
        viol_ = thermo_consis_fxn.getOutput('o0')
        return len(np.argwhere(viol_ < -tol)) == 0
    
    def CorrThermoConsis(self, dE_start, fix_Ea=[], fix_BE=[], print_screen=False):
        if not print_screen:
            old_stdout = sys.stdout
            sys.stdout = tempfile.TemporaryFile()

        Pnlp = self._Pnlp
        ini_p = np.hstack([dE_start])
        dev = Pnlp - ini_p

        object_fxn = cas.mul(dev.T, dev)
        nlp = dict(f=object_fxn, x=Pnlp, g=self._thermo_constraint_expression)

        nlpopts = dict()
        nlpopts['max_iter'] = 500
        nlpopts['tol'] = 1e-8
        nlpopts['acceptable_tol'] = 1e-8
        nlpopts['jac_d_constant'] = 'yes'
        nlpopts['expect_infeasible_problem'] = 'yes'
        nlpopts['hessian_approximation'] = 'exact'

        solver = cas.NlpSolver('solver', 'ipopt', nlp, nlpopts)

        # Bounds and initial guess
        lbP = -np.inf * np.ones(self._Np)
        ubP = np.inf * np.ones(self._Np)

        for i in range(len(fix_Ea)):
            if check_index(fix_Ea[i], self.dEa_index):
                lbP[i] = dE_start[find_index(fix_Ea[i], self.dEa_index)]
                ubP[i] = dE_start[find_index(fix_Ea[i], self.dEa_index)]
        for i in range(len(fix_BE)):
            if check_index(fix_BE[i], self.dBE_index):
                lbP[i + self._NEa] = dE_start[find_index(fix_BE[i], self.dBE_index) + len(self.dEa_index)]
                ubP[i + self._NEa] = dE_start[find_index(fix_BE[i], self.dBE_index) + len(self.dEa_index)]

        lbG = 0 * np.ones(2 * self.nrxn)
        ubG = np.inf * np.ones(2 * self.nrxn)

        solution = solver(x0=ini_p, lbg=lbG, ubg=ubG, lbx=lbP, ubx=ubP)
        dE_corr = solution['x'].full().T[0].tolist()

        if not print_screen:
            sys.stdout = old_stdout
        return(dE_corr)

    @property
    def Pnlp(self):
        return self._Pnlp


class KineticModel(SimpleKinetic):
    '''
    core kinetic model for forward simulation and parameter estimation
    '''
    def __init__(self, specieslist, reactionlist, dEa_index, dBE_index):
        SimpleKinetic.__init__(self, specieslist, reactionlist, dEa_index, dBE_index)
        self._x = None
        self._z = None
        self._p = None
        self._xdot = None
        self._zeq = None
        self._dae_ = None
        # Optimization
        self._evidence_ = None
        self._prior_ = None

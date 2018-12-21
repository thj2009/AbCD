# Define the class SPECIES, REACTION, REACTION NETWORK
from __future__ import print_function

import sys
sys.path.append(r'D:\Anaconda2\Lib\site-packages\casadi_2.4.2')
#sys.path.append('/home/hut216/CASADI/casadi_2.4.2')

import casadi as cas
import numpy as np
import scipy

# Physical Constant
c = 299792458           # Speed of light m s-1
kb = 1.38064852e-23     # Boltzmann Constant J K-1
h = 6.62607004e-34      # Blankc Constant   J s
NA = 6.02214086e23      # Avogadro Constant mol-1
R_g = 8.3144598         # Ideal Gas Constant
# Converter
hatree2kJmol = 2625.5   # Hatree to KJ.mol
hatree2J = 4.35974e-18  # Hatree to J
half_hc = 1/2. * h * c * 100 * NA/1000 # kJ cm /mol

class ActiveSite:
    def __init__(self, metal=''):
        self.metal = metal     
        self.struct = ''                # crystal structure
        self.mass = 0                   # mass of single atom, unit: g/mol
        self.facet = ''                 # exposed facet
        self.LatticeConstant = {}       # lattice constant of crystal structure
    
    def __repr__(self):
        return self.metal + '(' + self.facet + ')'
        
class Species:
    def __init__(self, name='', phase='', site=''):
        self.name = name
        self.phase = phase              # phase 'gas' or 'surface'
        self.site = None                # active site
        self.denticity = 0              # Number of site the surface species occupied
        self.element = {}
        self.BindInfo = None            # occupied site geometry configuration and binding energy
        self.mass = 0                   # mass of speceis unit: g/mol
        self.thermo = {'type':'', 'data': {}}
        self.vibfreq = []               # vibrational frequency
        self.BindEnergy = 0             # kJ/mol
        self.CorrectEnergy = 0          # Correction on Binding Energy 
        
    def __repr__(self):
        return self.name
    
    def readelement(self):
        spe = str(self.name)[0:-1]
        element = {'C': 0, 'H': 0, 'O': 0}
        if spe == 'Inert':
            return element
        for idx, ss in enumerate(spe):
            if not is_number(str(ss)):
                stoi = 1
                if not idx == len(spe) - 1:
                    if is_number(spe[idx+1]):
                        stoi = int(spe[idx+1])
                element[ss] += stoi
        return element
    

    def Enthalpy(self, Tem=273.15, correct=False):
        '''
        Calculate Enthapy or heat of formation at specific temperature
        Unit: kJ/mol
        '''
        H = 0
        if self.thermo['type'] == 'Shomate':
            shomate = self.thermo['data']
            tt = Tem/1000.

            if self.name != ''and bool(shomate):
                H = shomate['A']*tt + shomate['B']*tt**2/2. + shomate['C']*tt**3/3. \
                    + shomate['D']*tt**4/4. - shomate['E']/tt + shomate['F'] \
                    - shomate['H'] + shomate['dH']
            if correct:
                H += self.CorrectEnergy
                
        elif self.thermo['type'] == 'NASApoly':
            nasapoly = self.thermo['data']
            if self.name != '' and bool(nasapoly):
                H = R_g*Tem/1000 * (nasapoly['a1'] + nasapoly['a2']*Tem/2. + nasapoly['a3']*Tem**2/3. + nasapoly['a4']*Tem**3/4. + nasapoly['a5']*Tem**4/5. + nasapoly['a6']/Tem)
            if correct:
                H += self.CorrectEnergy
        return H
    
    def Entropy(self, Tem=273.15):
        S = 0
        if self.thermo['type'] == 'Shomate':
            shomate = self.thermo['data']
            tt = Tem/1000.
            if self.name != ''and bool(shomate):
                S = shomate['A']*np.log(tt) + shomate['B']*tt + shomate['C']*tt**2/2\
                    + shomate['D']*tt**3/3 - shomate['E']/(2*tt**2) + shomate['G']
        elif self.thermo['type'] == 'NASApoly':
            nasapoly = self.thermo['data']
            if self.name != ''and bool(nasapoly):
                S = R_g * (nasapoly['a1'] * np.log(Tem) + nasapoly['a2']*Tem + nasapoly['a3']*Tem**2/2. + nasapoly['a4']*Tem**3/3. + nasapoly['a5']*Tem**4/4. + nasapoly['a7'])            
        return S
    
    def Cp(self, Tem=273.15):
        Cp = 0
        if self.thermo['type'] == 'Shomate':
            shomate = self.thermo['data']
            tt = Tem/1000.
            if self.name != ''and bool(shomate):
                Cp = shomate['A'] + shomate['B']*tt + shomate['C']*tt**2 + shomate['D']*tt**3 - shomate['E']/(tt**2)
        elif self.thermo['type'] == 'NASApoly':
            nasapoly = self.thermo['data']
            if self.name != ''and bool(nasapoly):
                Cp = R_g * (nasapoly['a1'] + nasapoly['a2']*Tem + nasapoly['a3']*Tem**2 + nasapoly['a4']*Tem**3 + nasapoly['a5']*Tem**4)            
        return Cp
        
    def ZPEC(self):
        '''
        Calculate Zero point energy correction
        '''
        ZPE = 0
        for v in self.vibfreq:
            nablda = 1/v/100
            ZPE += 1/2.*h*c/nablda*NA/1000
        return ZPE

    def Entropy_vib(self, Tem=273.15):
        '''
        Calculate Vibrational Entropy
        '''
        Svib = 0
        for v in self.vibfreq:
            nablda = 0.01/v             # wavelength m
            x = h*c/nablda/(kb*Tem)
            Svib += x/(np.exp(x)-1) - np.log(1-np.exp(-x))
        Svib *= R_g                     # J mol-1 K-1
        return Svib
    
    def Entropy_2D(self, Tem=273.15):
        '''
        Calculate 2D gas Vibrational Entropy
        '''
        S2D = 0
        if self.phase == 'surface' and self.mass != 0:
            aa = self.site.LatticeConstant['a']
            ll = aa * np.sqrt(2)/2
            Area = np.sqrt(3)/4 * ll **2 * 2
            Strans = R_g * (np.log((2*np.pi*self.mass/1000/NA*kb*Tem)/h**2) + np.log(Area) + 2)
            Svib = 0
            for v in self.vibfreq[0:-2]:
                nablda = 0.01/float(v)      # wavelength m
                x = h*c/nablda/(kb*Tem)
                Svib += x/(np.exp(x)-1) - np.log(1-np.exp(-x))
            Svib *= R_g                     # J mol-1 K-1
            S2D = Strans + Svib
        else:
            S2D = 0
        return S2D
    
def ShomateFitting(Species):
    '''
    Fitting shomate parameter given species
    '''
    
class Reaction:
    '''
    Class Reaction
    represent the elementary reaction in surface science
    The reaction kinetics are calculated under the same temperature,
    Dynamic Temperature Add later
    '''
    def __init__(self, reactant = [], product = [], kinetic = None, name = None):
        # Note: Add name in the input file
        self.name = name            # name of reaction
        self.reactant = reactant    # list of reactant
        self.product = product      # list of product
        self.kinetic = kinetic      # Dictionary including kinetic type and parameter
        self.omega = 0              # Reaction Proximity
        self.deltaE = 0             # Delta E
        self.Ef = 0                 # Forward Reaction Barrier
        self.prefactor = 0
        self.CorrectEa = 0          # correct on Activation Barrier
#        self.CorrectTS = 0          # correct on transition state energy
        
    def __repr__(self):
        reactant = self.reactant
        product = self.product
        rxn = ''
        for i in range(len(reactant)):
            stoi = reactant[i][1]
            if stoi == 1:
                stoi_spe = str(reactant[i][0])
            else:
                stoi_spe = str(stoi) + str(reactant[i][0])
            if i == 0:
                rxn += stoi_spe
            else:
                rxn += '+' + stoi_spe
        rxn += '>>'
        for i in range(len(product)):
            stoi = product[i][1]
            if stoi == 1:
                stoi_spe = str(product[i][0])
            else:
                stoi_spe = str(stoi) + str(product[i][0])
            if i == 0:
                rxn += stoi_spe
            else:
                rxn += '+' + stoi_spe
        return rxn
    
    def Arrhenius(self, Tem=273.15, correct=False):
        Arr = {}
        if self.kinetic['type'] == 'Arr':
            Arr['A'] = self.kinetic['data']['A']
            Arr['Ea'] = self.kinetic['data']['Ea']
            Arr['n'] = self.kinetic['data']['n']
        elif self.kinetic['type'] == 'Shomate':
            H_R = self.IS_Enthalpy(Tem, correct)
            S_R = self.IS_Entropy(Tem)
            H_TS = self.TS_Enthalpy(Tem, correct)
            S_TS = self.TS_Entropy(Tem)
            Arr['A'] = kb/h*np.exp((S_TS - S_R)/R_g)
            Arr['Ea'] = H_TS - H_R
            Arr['n'] = 1
        return Arr
            
    def TS_Enthalpy(self, Tem=273.15, correct=False):
        if self.kinetic['type'] == 'Arr':
            ts = self.IS_Enthalpy(Tem, correct)
            ts += self.kinetic['data']['Ea']
            return ts
        elif self.kinetic['type'] == 'Shomate':
            shomate = self.kinetic['data']
            tt = Tem/1000.
            ts = 0
            if self.name != '':
                if bool(shomate):
                    ts = shomate['A']*tt + shomate['B']*tt**2/2. + shomate['C']*tt**3/3. \
                    + shomate['D']*tt**4/4. - shomate['E']/tt + shomate['F'] - shomate['H']\
                    + shomate['dH']
        if correct:
            corr_IS = self.IS_Enthalpy(Tem, 'True') - self.IS_Enthalpy(Tem, False)
            ts += self.CorrectEa + corr_IS
        return ts
    
    def TS_Entropy(self, Tem=273.15):
        '''
        Calculate the entropy of transition state using shomate parameter
        '''
        shomate = self.kinetic['data']
        tt = Tem/1000.
        S = 0
        if self.name != '':
            if bool(shomate):
                S = shomate['A']*np.log(tt) + shomate['B']*tt + shomate['C']*tt**2/2.\
                + shomate['D']*tt**3/3. - shomate['E']/(2*tt**2) + shomate['G']
        return S
    
    def IS_Enthalpy(self, Tem=273.15, correct=False):
        '''
        Calculate the enthalpy of Initial state using shomate parameter
        '''
        reactant = self.reactant
        H = 0
        for item in reactant:
            H += item[1]*item[0].Enthalpy(Tem, correct)
        return H
    
    def IS_Entropy(self, Tem=273.15):
        '''
        Calculate the enthalpy of Initial state using shomate parameter
        '''
        reactant = self.reactant
        S = 0
        for item in reactant:
            S += item[1]*item[0].Entropy(Tem)
        return S
    
    def FS_Enthalpy(self, Tem=273.15, correct=False):
        '''
        Calculate the entropy of Final state using shomate parameter
        '''
        product = self.product
        H = 0
        for item in product:
            H += item[1]*item[0].Enthalpy(Tem, correct)
        return H
    
    def FS_Entropy(self, Tem=273.15):
        '''
        Calculate the entropy of Final state using shomate parameter
        '''
        product = self.product
        S = 0
        for item in product:
            S += item[1]*item[0].Entropy(Tem)
        return S
    
    def dH_Enthalpy(self, Tem=273.15, correct=False):
        return self.FS_Enthalpy(Tem, correct) - self.IS_Enthalpy(Tem, correct)
    
    def dS_Entropy(self, Tem=273.15):
        return self.FS_Entropy(Tem) - self.IS_Entropy(Tem)
    
    def dG_Gibbs(self, Tem=273.15, correct=False):
        dH = self.dH_Enthalpy(Tem, correct)
        dS = self.dS_Entropy(Tem)
        dG = dH - Tem/1000*dS
        return dG
    
    def deltaH(self):
        '''
        Calculate the change on enthalpy given changing on species
        '''
        reactant = self.reactant
        H = 0
        for item in reactant:
            H -= item[1]*item[0].CorrectEnergy
        product = self.product
        for item in product:
            H += item[1]*item[0].CorrectEnergy
        return H
        
class ReactionNet:
    def __init__(self, specieslist = [], reactionlist = [], stoimat = []):
        self.specieslist = specieslist
        self.reactionlist = reactionlist
        self.Nspe = len(specieslist)
        self.Nrxn = len(reactionlist)
        self.stoimat = stoimat
        self.Ngas = self._Ngas()
        self.Nsurf = self._Nsurf()
        
    def __repr__(self):
        rxnnet = ''
        for rxn in self.reactionlist:
            rxnnet += str(rxn) + '\n'
        return rxnnet
    
    def _Ngas(self):
        Ngas = 0
        for spe in self.specieslist:
            if spe.phase == 'gaseous':
                Ngas += 1
        return Ngas
    
    def _Nsurf(self):
        Nsurf = 0
        for spe in self.specieslist:
            if spe.phase == 'surface':
                Nsurf += 1
        return Nsurf
        
    def StoiMat(self):
        stoimat = np.zeros([self.Nrxn, self.Nspe], dtype = 'int8')
        for i in range(self.Nrxn):
            # Deal with reactant
            reactant = self.reactionlist[i].reactant
            for rr in reactant:
                kk = get_index_species(str(rr[0]), self.specieslist)
                stoimat[i, kk] = -rr[1]
            prod = self.reactionlist[i].product
            for rr in prod:
                kk = get_index_species(str(rr[0]), self.specieslist)
                stoimat[i, kk] = rr[1]
        self.stoimat = stoimat
    
    def printscreen(self, file=''):
        if file != '':
            fout = open(file, 'w')
            old_stdout = sys.stdout
            sys.stdout = fout
        print('Reaction Network Summary')
        print(10 * self.Nspe * '=')
        print('Number of Species   : %d' %self.Nspe)
        print('Number of Reactions : %d' %self.Nrxn)
        print(10 * self.Nspe * '=')
        print('{:25s}'.format(' '), end = '')
        for i in range(self.Nspe):
            print('{:^10s}'.format(self.specieslist[i]), end = '')
        print('')
        for j in range(self.Nrxn):
            print('{:<25s}'.format(self.reactionlist[j]), end = '')
            for i in range(self.Nspe):
                print('{:^10d}'.format(self.stoimat[j, i]), end = '')
            print('')
        if file != '':
            fout.close()
            sys.stdout = old_stdout

class ReactionCondition:
    def __init__(self, name=''):
        self.name = name
        self.Temperature = 298.15       # unit Kelvin
        self.TotalPressure = 1          # total pressure, unit: atm
        self.TotalFlow = 0.00           # unit: s-1
        self.TotalSite = 0.00           # unit: total mole of activesite
        self.PartialPressure = {}        # Dictionary to store the partial pressure of gas phase species
        self.TurnOverFrequency = {}      # Dictionary to store outlet turnover frequency
        self.SimulationTime = 0
        self.SimulatedTOF = {}
        
    def __repr__(self):
        return self.name
    
    def printscreen(self, title=False):
        if title:
            print('{0:10s} {1:^15s} {2:^15s}'.format('name', 'Tem(K)', 'Ptot(atm)'), end = '')
            for key in self.PartialPressure.keys():
                print(' {0:^15s}'.format(key+'(atm)'), end = '')
            for key in self.TurnOverFrequency.keys():
                print(' {0:^15s}'.format(key+'(s-1)'), end = '')
            print()
        print('{0:10s} {1:^15.2f} {2:^15.2f}'.format(self.name, self.Temperature, self.TotalPressure), end = '')
        for val in self.PartialPressure.values():
            print(' {0:^15.2f}'.format(val), end = '')
        for val in self.TurnOverFrequency.values():
            print(' {0:^15.2e}'.format(val), end = '')
        print()
        
        
class CSTR_model:
    def __init__(self, reaction_net, dEa_index = [], dBE_index = []):
        self.reaction_net = reaction_net
        self.reactor = {}
        self.dEa_index = dEa_index
        self._dEa = self._dEa_f()
        self.dBE_index = dBE_index
        self._dBE = self._dBE_f()
        self._Np = self._Np_f()
        self._NEa = self._NEa_f()
        self._NBE = self._NBE_f()
        self._Pnlp = self._Pnlp_f()
        self._Tem = {}
        self._x = {}
        self._p = {}
        self._DAE = {}        
        self._thermo_constraint_expression = {}
        self._reaction_energy_expression = {}
        self._rate_expression = {}
        self._K_expression = {}
        self.TOF = 0
        self.partialpress = []
        self.coverage = []
        self.EqCons_RCons = {}
        self.Rate = {}                          # Store reaction rate information
        self.Energy = {}                        # Store reaction energy information
        self._dH_expression = {}
        self._Jacobian = {}
        
    def __repr__(self):
        return 'CSTRmodel'
    
    def _dEa_f(self):
        return cas.SX.sym('dEa', len(self.dEa_index))
    
    def _NEa_f(self):
        return len(self.dEa_index)
    
    def _NBE_f(self):
        return len(self.dBE_index)

    def _dBE_f(self):
        return cas.SX.sym('dBE', len(self.dBE_index))
        
    def _Np_f(self):
        return len(self.dBE_index) + len(self.dEa_index)
    
    def _Pnlp_f(self):
        #return cas.SX.sym('Pnlp', self._Np)
        return cas.MX.sym('Pnlp', self._Np)
    
    def Build_model(self):
        '''
        Create Ordinary Differential Equation for dynamic CSTR model
        And build the expression for energy, reaction kinetic, equilirium constant and 
        store in the class to call
        '''
        dEa = self._dEa
        dBE = self._dBE
        Pnlp = self._Pnlp
        
        stoimat = self.reaction_net.stoimat
        Nspe = self.reaction_net.Nspe
        Nrxn = self.reaction_net.Nrxn
        Ngas = self.reaction_net.Ngas
        Nsurf = self.reaction_net.Nsurf
        spelist = self.reaction_net.specieslist
        rxnlist = self.reaction_net.reactionlist


        # Define state variable
        flow = cas.SX.sym('flow', Ngas)                     # s-1
        cover = cas.SX.sym('cover', Nsurf)                  # fraction
        d_flow = cas.SX.sym('d_flow', Ngas)                 # s-2
        d_cover = cas.SX.sym('d_cover', Nsurf)              # fraction/s
        
        Ppress_inlet = cas.SX.sym('Ppress_inlet', Ngas)       # Inlet Partial Pressure of all gas phase species
        Ptot = 0                                              # Total Pressure 
        for i in range(Ngas):
            Ptot += Ppress_inlet[i]
        Ppress = cas.SX.sym('Ppress', Ngas)                   # Partial Pressure of gas species
        Tem = cas.SX.sym('Tem', 1)                            # Reaction Tempertature
        Flowtot = cas.SX.sym('Flowtot', 1)                    # total flowrate (s-1)
        for i in range(Ngas):
            Ppress[i] = flow[i]/Flowtot * Ptot
        

        #Calaulate reaction enthalpy
        enthal = []             # Store the enthalpy of each species to build DAE
        enthal_cons = []        # Store the enthalpy of each species to build constraint
        entro = []
        e0 = []                 # Store the enthalpy of each species for function calculation
        denthal = []
        for j in range(Nspe):
            HH = spelist[j].Enthalpy(Tem)
            SS = spelist[j].Entropy(Tem)
            enthal.append(HH) 
            enthal_cons.append(spelist[j].Enthalpy(298.15))
            entro.append(SS)
            e0.append(HH)
            denthal.append(0)
        
        for j in range(len(self.dBE_index)):          # surface species add deviation variable
            enthal[self.dBE_index[j]] += dBE[j]
            enthal_cons[self.dBE_index[j]] += Pnlp[len(self.dEa_index)+j]
            denthal[self.dBE_index[j]] += Pnlp[len(self.dEa_index)+j]
            
        # calculate reaction kinetic constants (forward and reverse) of each reactions
        kf = cas.SX.sym('kf', Nrxn)
        kr = cas.SX.sym('kr', Nrxn)
        Keq = cas.SX.sym('Keq', Nrxn)
        enthal_react_cons = []          # Store the Enthalpy of each reaction
        Ea_cons =[]                     # Store the Activation Barriar of each reaction
        enthal_react = []; entro_react = []
        react0 = []; Ea_00 = []         # Store the reaction for function calculation
        denthal_react = []
        for i in range(Nrxn):
            kine_data = rxnlist[i].Arrhenius(Tem)
            Hreact = 0; dHreact = 0
            Sreact = 0
            Hreact_cons = 0; H0 = 0
            Arr = kine_data['A']
            Ea0 = kine_data['Ea']
            Ea298 = rxnlist[i].Arrhenius(298.15)['Ea']
            n = kine_data['n']
            if check_index(i, self.dEa_index):
                ind = find_index(i, self.dEa_index)
                Ea = Ea0 + dEa[ind]
                Ea_cons.append(Ea298 + Pnlp[ind])
            else:
                Ea = Ea0
                Ea_cons.append(Ea298)
            Ea_00.append(Ea)
            # k = A * T^n * exp(-Ea/(RT))
            kf[i] = Arr * cas.exp(-Ea*1000/(R_g*Tem)) * (Tem**n)
            for j in range(Nspe):
                Hreact += stoimat[i][j]*enthal[j]
                dHreact += stoimat[i][j]*denthal[j]
                Hreact_cons += stoimat[i][j]*enthal_cons[j]
                Sreact += stoimat[i][j]*entro[j]
                H0 += stoimat[i][j]*e0[j]
    
            enthal_react.append(Hreact)
            denthal_react.append(dHreact)
            enthal_react_cons.append(Hreact_cons)
            entro_react.append(Sreact)
            react0.append(H0)
    
            Keq[i] = cas.exp(-Hreact*1000/(R_g*Tem)) * cas.exp(Sreact/R_g)
            kr[i] = kf[i]/Keq[i]
        
        rate = cas.SX.sym('rate', Nrxn)         # s-1
        rfor = cas.SX.sym('rfor', Nrxn)         # s-1
        rrev = cas.SX.sym('rrev', Nrxn)         # s-1
        Qeq = cas.SX.sym('Qeq', Nrxn)
        for i in range(Nrxn):
            rfor[i]=kf[i]; rrev[i]=kr[i]; Qeq[i] = 1
            for j in range(Nspe):
                if stoimat[i][j] < 0:
                    if j < Ngas:
                        rfor[i] *= Ppress[j]**(-stoimat[i][j])
                        Qeq[i] /= Ppress[j]**(-stoimat[i][j])
                    else:
                        rfor[i] *= (cover[j-Ngas])**(-stoimat[i][j])
                        Qeq[i] /= (cover[j-Ngas])**(-stoimat[i][j])
                elif stoimat[i][j] > 0:
                    if j < Ngas:
                        rrev[i] *= Ppress[j]**(stoimat[i][j])
                        Qeq[i] *= Ppress[j]**(stoimat[i][j])
                    else:
                        rrev[i] *= (cover[j-Ngas])**(stoimat[i][j])
                        Qeq[i] *= (cover[j-Ngas])**(stoimat[i][j])
            rate[i] = rfor[i] - rrev[i]
        
        tau = 1
        for j in range(Ngas):
            d_flow[j] = Ppress_inlet[j] / Ptot * Flowtot - flow[j]
            for i in range(Nrxn):
                d_flow[j] += stoimat[i][j]*rate[i]
            d_flow[j] /= tau
                       
        for j in range(Nsurf):
            d_cover[j] = 0
            for i in range(Nrxn):
                d_cover[j] += stoimat[i][j+Ngas]*rate[i]
        
        x = cas.vertcat([flow, cover])
        p = cas.vertcat([dEa, dBE, Ppress_inlet, Tem, Flowtot])
        x_dot = cas.vertcat([d_flow, d_cover])
        
        dae = dict(x=x, p=p, ode=x_dot)
        
        #create constraint on Parameter
        ineq2 = cas.vertcat(Ea_cons)-cas.vertcat(enthal_react_cons)
        thermal_consis_ineq = cas.vertcat([cas.vertcat(Ea_cons), ineq2]);                   
        self._x = x
        self._p = p
        self._DAE = dae
        self._Tem = Tem
        self._thermo_constraint_expression = thermal_consis_ineq
        self._reaction_energy_expression = {'enthalpy': cas.vertcat(enthal_react), 'activation': cas.vertcat(Ea_00)}
        self._rate_expression = {'rnet':rate, 'rfor': rfor, 'rrev': rrev}
        self._K_expression = {'Keq': Keq, 'Qeq': Qeq, 'kf': kf, 'kr': kr}
        self._dH_expression = denthal_react
        
    def ForwardCal(self, dEa_start, dBE_start, Condition, Fullcal = False, abstol=1e-12, reltol=1e-10, DRC=False):
        Ngas = self.reaction_net.Ngas
        TotalPressure = Condition.TotalPressure 
        TotalFlow = Condition.TotalFlow
        Tem = Condition.Temperature
        tf = Condition.SimulationTime
        opts = {}
        opts['tf'] = tf       # Simulation time
        opts['abstol'] = abstol
        opts['reltol'] = reltol
        opts['linear_multistep_method'] = 'bdf'
        opts['disable_internal_warnings'] = True
        opts['max_num_steps'] = 1e5
    
        Fint = cas.Integrator('Fint', 'cvodes', self._DAE, opts)
        
        x0 = [0] * (self.reaction_net.Nspe - 1) + [1]    
        
        # Partial Pressure
        Pinlet = np.zeros(Ngas)
        for idx, spe in enumerate(self.reaction_net.specieslist):
            if spe.phase == 'gaseous':
                Pinlet[idx] = Condition.PartialPressure[str(spe)]
        P_dae = np.hstack([dEa_start, dBE_start, Pinlet, Tem, TotalFlow])
        F_sim = Fint(x0=x0, p=P_dae)
        
        tor = {}
        for idx, spe in enumerate(self.reaction_net.specieslist):
            if spe.phase == 'gaseous':
                tor[str(spe)] = float(Pinlet[idx]/TotalPressure * TotalFlow - F_sim['xf'][idx])
#        print(Pinlet)
#        print(TotalPressure, TotalFlow)
#        print(F_sim['xf'])
        
        self.TOF = tor
        self.partialpress = (F_sim['xf'][:Ngas]/TotalFlow * TotalPressure).full().T[0].tolist()
        self.coverage = F_sim['xf'][Ngas:].full().T[0].tolist()
        
        # May Create Individual function to call this
        # Evaluate Reaction Rate automatically save to Rate attribute
        x = self._x
        p = self._p
        RateFxn = cas.SXFunction('RateFxn', [x, p], [self._rate_expression['rnet'], self._rate_expression['rfor'], self._rate_expression['rrev']])
        RateFxn.setInput(F_sim['xf'], 'i0')
        RateFxn.setInput(P_dae, 'i1')
        RateFxn.evaluate()
        self.Rate['rnet'] = RateFxn.getOutput('o0').full().T[0].tolist()
        self.Rate['rfor'] = RateFxn.getOutput('o1').full().T[0].tolist()
        self.Rate['rrev'] = RateFxn.getOutput('o2').full().T[0].tolist()
        
        # Evaluate Reaction Energy
        EnFxn = cas.SXFunction('EnFxn', [x, p], [self._reaction_energy_expression['activation'], self._reaction_energy_expression['enthalpy']])
        EnFxn.setInput(F_sim['xf'], 'i0')
        EnFxn.setInput(P_dae, 'i1')
        EnFxn.evaluate()
        self.Energy['activation'] = EnFxn.getOutput('o0').full().T[0].tolist()
        self.Energy['enthalpy'] = EnFxn.getOutput('o1').full().T[0].tolist()
        
        # Evaluate Equilibrium Constant and Rate Constant
        KFxn = cas.SXFunction('KFxn', [x, p], [self._K_expression['Keq'], self._K_expression['Qeq'], self._K_expression['kf'], self._K_expression['kr']])
        KFxn.setInput(F_sim['xf'], 'i0')
        KFxn.setInput(P_dae, 'i1')
        KFxn.evaluate()
        self.EqCons_RCons['Keq'] = KFxn.getOutput('o0').full().T[0].tolist()
        self.EqCons_RCons['Qeq'] = KFxn.getOutput('o1').full().T[0].tolist()
        self.EqCons_RCons['kf'] = KFxn.getOutput('o2').full().T[0].tolist()
        self.EqCons_RCons['kr'] = KFxn.getOutput('o3').full().T[0].tolist()
        
        Condition.SimulatedTOF = tor
        
        if DRC:
            opts = {}
            opts['tf'] = tf      # Simulation time
            opts['abstol'] = 1e-10
            opts['reltol'] = 1e-8
            opts['abstolB'] = 1e-6
            opts['reltolB'] = 1e-4
            opts['linear_multistep_method'] = 'bdf'
            opts['disable_internal_warnings'] = True
            opts['max_num_steps'] = 1e5
            Fint = cas.Integrator('Fint', 'cvodes', self._DAE, opts)
            
            
            Pnlp = self._Pnlp
            P_dae = cas.vertcat([Pnlp, Pinlet, Tem, TotalFlow])
            #print(P_dae)
            F_sim = Fint(x0=x0, p=P_dae)
            
            XRC = {}
            for idx, spe in enumerate(self.reaction_net.specieslist):
                if spe.phase == 'gaseous':
                    # TOR
                    tor_DRC = Pinlet[idx]/TotalPressure * TotalFlow - F_sim['xf'][idx]
                    # Jacobian
                    Jac = cas.jacobian(tor_DRC, Pnlp)
                    F_jac = cas.MXFunction('F_jac', [P_], [Jac, tor_DRC])
                    PP = np.hstack([dEa_start, dBE_start])
                    F_jac.setInput(PP, 'i0')
                    F_jac.evaluate()
                    JJ = F_jac.getOutput('o0').full()[0]
                    #print(JJ)
                    Tor = F_jac.getOutput('o1').full()
                    #print(Tor)
                    XRC[str(spe)] = -R_g*Tem/Tor*JJ
            self.DRC = XRC
        return tor
        
        
    def Single_Opt(self, dEa_start, dBE_start, ConditionList, L2=1e-5,
                   Thermo=True, Printscreen=False, Report=None,
                   tol=1e-2, maxiter=5000):
        if not Printscreen:
            import tempfile
            old_stdout = sys.stdout
            sys.stdout = tempfile.TemporaryFile()
        elif Report:
            old_stdout = sys.stdout
            sys.stdout = open(Report, 'w')
            
        
        Pnlp = self._Pnlp
        
        Ngas = self.reaction_net.Ngas
        opts = optimize_option()
        Fint = cas.Integrator('Fint', 'cvodes', self._DAE, opts)
        
        x0 = [0] * (self.reaction_net.Nspe - 1) + [1]    
        
        obj = cas.mul(Pnlp.T, Pnlp)*L2
        for Condition in ConditionList:
            TotalPressure = Condition.TotalPressure 
            TotalFlow = Condition.TotalFlow
            Tem = Condition.Temperature

            # Partial Pressure
            Pinlet = np.zeros(Ngas)
            for idx, spe in enumerate(self.reaction_net.specieslist):
                if spe.phase == 'gaseous':
                    Pinlet[idx] = Condition.PartialPressure[str(spe)]
            P_dae = cas.vertcat([Pnlp, Pinlet, Tem, TotalFlow])
            F_sim = Fint(x0=x0, p=P_dae)
            for idx, spe in enumerate(self.reaction_net.specieslist):
                if spe.phase == 'gaseous':
                    # TOR
                    tor = Pinlet[idx]/TotalPressure * TotalFlow - F_sim['xf'][idx]
                    if str(spe) in Condition.TurnOverFrequency.keys():
                        err = tor - Condition.TurnOverFrequency[str(spe)]
                        obj += err * err
                        #print(err)


        nlp = dict(f=obj, x=Pnlp, g=self._thermo_constraint_expression)
        
        nlpopts = {}
        nlpopts['max_iter'] = maxiter
        nlpopts['tol'] = tol
        nlpopts['acceptable_tol'] = 1e-3
        nlpopts['jac_d_constant'] = 'yes'
        nlpopts['expect_infeasible_problem'] = 'yes'
        nlpopts['hessian_approximation'] = 'limited-memory'
        nlpopts['gather_stats'] = True
        nlpopts['print_level'] = 5
        solver = cas.NlpSolver('solver', 'ipopt', nlp, nlpopts)
        
        # Bounds and initial guess
        lbP = -40 * np.ones(self._Np)
        ubP = 40 * np.ones(self._Np)
        
        lbG = 0 * np.ones(2 * self.reaction_net.Nrxn)
        ubG = np.inf * np.ones(2 * self.reaction_net.Nrxn)
        
        x0 = np.hstack([dEa_start, dBE_start])
        solution = solver(x0=x0, lbx=lbP, ubx=ubP, lbg=lbG, ubg=ubG)   
        opt_sol = solution['x'].full().T[0].tolist()
        obj = solution['f'].full()[0][0]
        
        print('===============================================================')
        print('Starting Point:')
        print(dEa_start)
        print(dBE_start)
        print('Parameter:')
        print(opt_sol)
        print('Objective:')
        print(obj)
        
#        
#        print('**************************')
#        # Evaluate Objective function
#        objfxn = cas.MXFunction('objfxn', [Pnlp], [obj, Tor_s])
#        objfxn.setInput(opt_sol, 'i0')
#        objfxn.evaluate()
#        print(objfxn.getOutput('o0'))
#        print(objfxn.getOutput('o1'))
        if not Printscreen:
            sys.stdout = old_stdout
        elif Report:
            sys.stdout.close()
            sys.stdout = old_stdout

        return(opt_sol, obj)
        
        
    def CheckThermoConsis(self, dEa_start, dBE_start, tol = 1e-8): 
        '''
        True: satify the thermodynamic consistency
        False: not satisfy
        '''
        Pnlp = self._Pnlp
        #Tem = self._Tem
        ThermoConsisFxn = cas.MXFunction('ThermoConsisFxn', [Pnlp], [self._thermo_constraint_expression])
        Ini_p = np.hstack([dEa_start, dBE_start])
        #print(Ini_p.shape)
        #print(Pnlp.shape)
        #print(self._thermo_constraint_expression.shape)
        ThermoConsisFxn.setInput(Ini_p, 'i0')
        ThermoConsisFxn.evaluate()
        viol_ = ThermoConsisFxn.getOutput('o0')
#        print(viol_)
#        print(len(np.argwhere(viol_ < -tol)))
        return len(np.argwhere(viol_ < -tol)) == 0
            
        
    def CorrThermoConsis(self, dEa_start, dBE_start, fix_Ea = [], fix_BE = [], Corr = True, Printscreen = False):
        import tempfile
        if not Printscreen:
            old_stdout = sys.stdout
            sys.stdout = tempfile.TemporaryFile()
        
        Pnlp = self._Pnlp
        #Tem = self._Tem
        Ini_p = np.hstack([dEa_start, dBE_start])
        dev = Pnlp - Ini_p
        
        object_fxn = cas.mul(dev.T, dev)
        nlp = dict(f=object_fxn, x = Pnlp, g=self._thermo_constraint_expression)
        
        nlpopts = dict()   
        nlpopts['max_iter'] = 500
        nlpopts['tol'] = 1e-8
        nlpopts['acceptable_tol'] = 1e-8
        nlpopts['jac_d_constant'] = 'yes'
        nlpopts['expect_infeasible_problem'] = 'yes'
        nlpopts['hessian_approximation'] = 'exact'
        nlpopts['print_level'] = 5
        solver = cas.NlpSolver('solver', 'ipopt', nlp, nlpopts)
        
    #    # Bounds and initial guess
        lbP = -np.inf * np.ones(self._Np)
        ubP = np.inf * np.ones(self._Np)
        
        for i in range(len(fix_Ea)):
            if check_index(fix_Ea[i], self.dEa_index):
                lbP[i] = dEa_start[find_index(fix_Ea[i], self.dEa_index)] 
                ubP[i] = dEa_start[find_index(fix_Ea[i], self.dEa_index)]
        for i in range(len(fix_BE)):
            if check_index(fix_BE[i], self.dBE_index):
                lbP[i + self._NEa] = dBE_start[find_index(fix_BE[i], self.dBE_index)]
                ubP[i + self._NEa] = dBE_start[find_index(fix_BE[i], self.dBE_index)]

        lbG = 0 * np.ones(2 * self.reaction_net.Nrxn)
        ubG = np.inf * np.ones(2 * self.reaction_net.Nrxn)
        
        solution = solver(x0 = Ini_p, lbg = lbG, ubg = ubG, lbx = lbP, ubx = ubP)    
        xx = solution['x'].full().T[0].tolist()
        
        if not Printscreen:
            sys.stdout = old_stdout
        dEa_corr = xx[:self._NEa]; dBE_corr = xx[self._NEa:]
        
        return(dEa_corr, dBE_corr)


    def Single_Opt_prior(self, dEa_start, dBE_start, Full_Pinlet, Full_Measure,
                   bound = [-20, 20], Mean_BE = [], invCov_BE = [], Sig_BE=[], 
                   Mean_Ea=[], Sig_Ea=[], 
                   err_obs = 0.1, 
                   Thermo = True, Print_opt = False, Report=None,
                   tol = 1e-1, maxiter = 500):
        if not Print_opt:
            import tempfile
            old_stdout = sys.stdout
            sys.stdout = tempfile.TemporaryFile()
        elif Report:
            old_stdout = sys.stdout
            sys.stdout = open(Report, 'w')
        Pnlp = self._Pnlp
        reactor = self.reactor
        # Define Experiment Parameter
        inlet_flow = reactor['inletflow'] # m^3/s
        inlet_mole = reactor['inletmole']
        inlet_press = reactor['inletpress']
        Tem = reactor['temperature']  # Kelvin
        cata_mole = reactor['cata_mole']

        opts = optimize_option()
        Fint = cas.Integrator('Fint', 'cvodes', self._DAE, opts)
        
        x0 = [0] * (self.reaction_net.Nspe - 1) + [1] 
        Tor_s = []
        
        Nexp = len(Full_Measure)
        for k in range(Nexp):
            Pinlet = Full_Pinlet[:,k]; 
            P_dae = cas.vertcat([Pnlp,  Pinlet])
            F_sim = Fint(x0=x0, p=P_dae)
            tor= (Pinlet[1]/inlet_press*inlet_mole - F_sim['xf'][1])/cata_mole*100             
            Tor_s.append(tor)
            
        Tor_s = cas.vertcat(Tor_s)
        err = Tor_s - np.array(Full_Measure)
        
        # Deviation on Turn over rate
        obj_deviation = cas.mul(err.T, err)/(2*err_obs**2)
        
        # Deviation on Surface Energy
        ddBE = Pnlp[self._NEa:] - Mean_BE
        obj_BE = 1/2.*cas.mul(cas.mul(ddBE.T, invCov_BE), ddBE)
        
        # Deviation on Activation Barrier
        obj_Ea = 0
        for i in range(self._NEa):
            obj_Ea += ((Pnlp[i] - Mean_Ea[i])**2)/(2*Sig_Ea[i]**2)
            
        obj = obj_deviation + obj_BE + obj_Ea
        nlp = dict(f=obj, x=Pnlp, g=self._thermo_constraint_expression)
        
        nlpopts = {}
        nlpopts['max_iter'] = maxiter
        nlpopts['tol'] = tol
        nlpopts['acceptable_tol'] = 1e-3
        nlpopts['jac_d_constant'] = 'yes'
        nlpopts['expect_infeasible_problem'] = 'yes'
        nlpopts['hessian_approximation'] = 'limited-memory'
        nlpopts['gather_stats'] = True
        nlpopts['print_level'] = 5
        solver = cas.NlpSolver('solver', 'ipopt', nlp, nlpopts)
        
        # Bounds and initial guess
        lbP = bound[0] * np.ones(self._Np)
        ubP = bound[1] * np.ones(self._Np)
        lbG = 0 * np.ones(2 * self.reaction_net.Nrxn)
        ubG = np.inf * np.ones(2 * self.reaction_net.Nrxn)

        x0 = np.hstack([dEa_start, dBE_start])
        solution = solver(x0=x0, lbx=lbP, ubx=ubP, lbg=lbG, ubg=ubG)   
        opt_sol = solution['x'].full().T[0].tolist()
        obj = solution['f'].full()[0][0]
        
        print('===============================================================')
        print('Starting Point:')
        print(dEa_start)
        print(dBE_start)
        print('Parameter:')
        print(opt_sol)
        print('Objective:')
        print(obj)
        
#        print('**************************')
#        # Evaluate Objective function
#        objfxn = cas.MXFunction('objfxn', [Pnlp], [obj_deviation, obj_BE, obj_omega])
#        objfxn.setInput(opt_sol, 'i0')
#        objfxn.evaluate()
#        print(objfxn.getOutput('o0'))
#        print(objfxn.getOutput('o1'))
#        print(objfxn.getOutput('o2'))
#
#        # Evaluate Objective function
#        aafxn = cas.MXFunction('objfxn', [Pnlp], [ddomega])
#        aafxn.setInput(opt_sol, 'i0')
#        aafxn.evaluate()
#        print(aafxn.getOutput('o0'))

        if not Print_opt:
            sys.stdout = old_stdout
        elif Report:
            sys.stdout.close()
            sys.stdout = old_stdout
        return(opt_sol, obj)

    def Single_Opt_prior_omega(self, dEa_start, dBE_start, Full_Pinlet, Full_Measure,
                   bound = [-20, 20], Mean_BE = [], invCov_BE = [], Sig_BE=[],
                   err_obs = 0.1, Sig_Ea = [],
                   Thermo = True, Print_opt = False, Report=None,
                   tol = 1e-2, maxiter = 500):
        if not Print_opt:
            import tempfile
            old_stdout = sys.stdout
            sys.stdout = tempfile.TemporaryFile()
        elif Report:
            old_stdout = sys.stdout
            sys.stdout = open(Report, 'w')
        Pnlp = self._Pnlp
        reactor = self.reactor
        Rxn_list = self.reaction_net.reactionlist
        # Define Experiment Parameter
        inlet_flow = reactor['inletflow'] # m^3/s
        inlet_mole = reactor['inletmole']
        inlet_press = reactor['inletpress']
        Tem = reactor['temperature']  # Kelvin
        cata_mole = reactor['cata_mole']

        opts = optimize_option()
        Fint = cas.Integrator('Fint', 'cvodes', self._DAE, opts)
        
        x0 = [0] * (self.reaction_net.Nspe - 1) + [1] 
#        x0 = [1e-6]*12 + [1]
        Tor_s = []
        
        Nexp = len(Full_Measure)
        for k in range(Nexp):
            Pinlet = Full_Pinlet[:,k]; 
            P_dae = cas.vertcat([Pnlp,  Pinlet])
            F_sim = Fint(x0=x0, p=P_dae)
            tor= (Pinlet[1]/inlet_press*inlet_mole - F_sim['xf'][1])/cata_mole*100             
            Tor_s.append(tor)
            
        Tor_s = cas.vertcat(Tor_s)
        err = Tor_s - np.array(Full_Measure)
        # Deviation on Turn over rate
        obj_deviation = cas.mul(err.T, err)/(2*err_obs**2)
        
        # Deviation on Surface Energy
        ddBE = Pnlp[self._NEa:] - Mean_BE
        obj_BE = 1/2.*cas.mul(cas.mul(ddBE.T, invCov_BE), ddBE)
#        obj_BE = 0
#        for i in range(self._NBE):
#            obj_BE += ((Pnlp[self._NEa+i] - Mean_BE[i])**2)/(2*Sig_BE[i]**2)
#        print(obj_BE)
        # Deviation on Surface Energy
        obj_Ea = 0; 
        for i in range(self._NEa):
            rxn = Rxn_list[self.dEa_index[i]]
            dEa_exp = rxn.omega * self._dH_expression[self.dEa_index[i]]
            ddEa = Pnlp[i] - dEa_exp
            obj_Ea += (ddEa**2)/(2*Sig_Ea**2)

        # Sum all three 
        obj = (obj_deviation + obj_BE + obj_Ea) * (2*err_obs**2)
        nlp = dict(f=obj, x=Pnlp, g=self._thermo_constraint_expression)
        
        nlpopts = {}
        nlpopts['max_iter'] = maxiter
        nlpopts['tol'] = tol
        nlpopts['acceptable_tol'] = 1e-3
        nlpopts['jac_d_constant'] = 'yes'
        nlpopts['expect_infeasible_problem'] = 'yes'
        nlpopts['hessian_approximation'] = 'limited-memory'
        nlpopts['gather_stats'] = True
        nlpopts['print_level'] = 5
        solver = cas.NlpSolver('solver', 'ipopt', nlp, nlpopts)
        
        # Bounds and initial guess
        lbP = bound[0] * np.ones(self._Np)
        ubP = bound[1] * np.ones(self._Np)
        lbG = 0 * np.ones(2 * self.reaction_net.Nrxn)
        ubG = np.inf * np.ones(2 * self.reaction_net.Nrxn)
        x0 = np.hstack([dEa_start, dBE_start])
        
        solution = solver(x0=x0, lbx=lbP, ubx=ubP, lbg=lbG, ubg=ubG)   
        opt_sol = solution['x'].full().T[0].tolist()
        obj = solution['f'].full()[0][0]
        
        print('===============================================================')
        print('Starting Point:')
        print(dEa_start)
        print(dBE_start)
        print('Parameter:')
        print(opt_sol)
        print('Objective:')
        print(obj)

        if not Print_opt:
            sys.stdout = old_stdout
        elif Report:
            sys.stdout.close()
            sys.stdout = old_stdout
        return(opt_sol, obj)

    def SVD_Jacobian(self, Popt, Full_Pinlet, Full_Measure,
                   Mean_BE = [], Cov_BE = [],
                   Mean_Ea=[], Sig_Ea=[], 
                   err_obs = 0.1):
        '''
        Calculate the Jacobian of objective function at given point
        '''
        Pnlp = self._Pnlp
        reactor = self.reactor
        # Define Experiment Parameter
        inlet_flow = reactor['inletflow'] # m^3/s
        inlet_mole = reactor['inletmole']
        inlet_press = reactor['inletpress']
        Tem = reactor['temperature']  # Kelvin
        cata_mole = reactor['cata_mole']

        opts = optimize_option()
        Fint = cas.Integrator('Fint', 'cvodes', self._DAE, opts)
        
        x0 = [0] * (self.reaction_net.Nspe - 1) + [1] 
        Tor_s = []
        
        Nexp = len(Full_Measure)
        for k in range(Nexp):
            Pinlet = Full_Pinlet[:,k]; 
            P_dae = cas.vertcat([Pnlp,  Pinlet])
            F_sim = Fint(x0=x0, p=P_dae)
            tor= (Pinlet[1]/inlet_press*inlet_mole - F_sim['xf'][1])/cata_mole*100             
            Tor_s.append(tor)
#            dtor = cas.jacobian(tor, P_nlp)
        
        Tor_s = cas.vertcat(Tor_s)
        Jac_deviation = cas.jacobian(Tor_s/err_obs, Pnlp)
        CC = scipy.linalg.block_diag(np.diag(Sig_Ea**2), Cov_BE)
        H_prior = np.linalg.inv(CC)
        Jac_Fxn = cas.MXFunction('Jac_Fxn', [Pnlp], [Jac_deviation])
        Jac_Fxn.setInput(Popt, 'i0')
        Jac_Fxn.evaluate()
        JJ = Jac_Fxn.getOutput('o0')
        JJ = np.array(JJ)
        HH = (JJ.T).dot(JJ) + H_prior
        
#        # Exact Hessian 
#        err = Tor_s - np.array(Full_Measure)
#        # Deviation on Turn over rate
#        obj_deviation = cas.mul(err.T, err)/(2*err_obs**2)
#        # Deviation on Surface Energy
#        ddBE = Pnlp[self._NEa:] - Mean_BE
#        obj_BE = 1/2.*cas.mul(cas.mul(ddBE.T, np.linalg.inv(Cov_BE)), ddBE)
#        # Deviation on Activation Barrier
#        obj_Ea = 0
#        for i in range(self._NEa):
#            obj_Ea += ((Pnlp[i] - Mean_Ea[i])**2)/(2*Sig_Ea[i]**2)
#        obj = obj_deviation + obj_BE + obj_Ea
##        obj = obj_BE + obj_Ea
#        H_deviation = cas.hessian(obj, Pnlp)[0]
#        H_Fxn = cas.MXFunction('Jac_Fxn', [Pnlp], [H_deviation])
#        H_Fxn.setInput(Popt, 'i0')
#        H_Fxn.evaluate()
#        HH = H_Fxn.getOutput('o0')
        
        Eig, Vec = np.linalg.eig(HH)
        # HH = Vec * diag(Eig) * Vec.T
#        print(np.allclose(HH, (Vec).dot(np.diag(Eig)).dot(Vec.T)))
#        print(np.allclose(HH.T, HH))
        return(Eig, Vec)


    def SVD_Jacobian_omega(self, Popt, Full_Pinlet, Full_Measure,
                   Mean_BE = [], Cov_BE = [], Std_Ea=10,
                   err_obs = 0.05):
        '''
        Calculate the Jacobian of objective function at given point
        '''
        Pnlp = self._Pnlp
        reactor = self.reactor
        Rxn_list = self.reaction_net.reactionlist
        # Define Experiment Parameter
        inlet_flow = reactor['inletflow'] # m^3/s
        inlet_mole = reactor['inletmole']
        inlet_press = reactor['inletpress']
        Tem = reactor['temperature']  # Kelvin
        cata_mole = reactor['cata_mole']
        opts = optimize_option()
        Fint = cas.Integrator('Fint', 'cvodes', self._DAE, opts)
        
        x0 = [0] * (self.reaction_net.Nspe - 1) + [1] 
        Tor_s = []
        
        Nexp = len(Full_Measure)
        for k in range(Nexp):
            Pinlet = Full_Pinlet[:,k]; 
            P_dae = cas.vertcat([Pnlp,  Pinlet])
            F_sim = Fint(x0=x0, p=P_dae)
            tor= (Pinlet[1]/inlet_press*inlet_mole - F_sim['xf'][1])/cata_mole*100             
            Tor_s.append(tor)
        Tor_s = cas.vertcat(Tor_s)
        
        # Jacobian on TOR Deviation
        Jac_tor_deviation = cas.jacobian(Tor_s/err_obs, Pnlp)
        
        # Prior on Surface Energy
        HBE_prior = np.linalg.inv(Cov_BE);
        H_prior = scipy.linalg.block_diag(np.zeros([self._NEa, self._NEa]), HBE_prior)
        # Deviation on Surface Energy
        ddEa = []
        for i in range(self._NEa):
            rxn = Rxn_list[self.dEa_index[i]]
            dEa_exp = rxn.omega * self._dH_expression[self.dEa_index[i]]
            ddEa.append(Pnlp[i] - dEa_exp)
        ddEa = cas.vertcat(ddEa)
        Jac_Ea_deviation = cas.jacobian(ddEa/Std_Ea, Pnlp)
        # Jacobian on TOF
        Jac_tor_fxn = cas.MXFunction('Jac_tor_fxn', [Pnlp], [Jac_tor_deviation])
        Jac_tor_fxn.setInput(Popt, 'i0')
        Jac_tor_fxn.evaluate()
        JJ_tor = Jac_tor_fxn.getOutput('o0')
        JJ_tor = np.array(JJ_tor)
        # Jacobian on dEa
        Jac_Ea_fxn = cas.MXFunction('Jac_Ea_fxn', [Pnlp], [Jac_Ea_deviation])
        Jac_Ea_fxn.setInput(Popt, 'i0')
        Jac_Ea_fxn.evaluate()
        JJ_Ea = Jac_Ea_fxn.getOutput('o0')
        JJ_Ea = np.array(JJ_Ea)
        

        HH = (JJ_tor.T).dot(JJ_tor) + H_prior + (JJ_Ea.T).dot(JJ_Ea)
        Eig, Vec = np.linalg.eigh(HH)
        

        return(Eig, Vec)




def check_index(index, checklist):
    for i in checklist:
        if i == index:
            return True
    return False

def find_index(index, checklist):
    for i in range(len(checklist)):
        if index == checklist[i]:
            return i

def get_index_species(spe, spelist):
    for i in range(len(spelist)):
            if spe == str(spelist[i]):
                return i
    raise NameError(spe + ' is not in the species list')


def get_index_reactions(rxnname, rxnlist):
    for i in range(len(rxnlist)):
        if rxnname == rxnlist[i].name:
            return i
    raise NameError(rxnname+ ' is not in the reaction list')

def get_index_condition(condiname, condilist):
    for i in range(len(condilist)):
        if condiname == condilist[i].name:
            return i
    raise NameError(condiname+ ' is not in the condition list')
    
def optimize_option():
    '''
    Options pass to IPOPT
    '''
    tf = 5000
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
    opts['abstol'] = 1e-10
    opts['reltol'] = 1e-8
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

def AssignCorrection(specieslist, reactionlist, dBE_index=[], dEa_index=[], Change=[]):
    sol_Ea = Change[0:len(dEa_index)]
    sol_BE = Change[len(dEa_index):]

    for idx, key in enumerate(dBE_index):
        spe = specieslist[key]
        spe.CorrectEnergy = sol_BE[idx]
    for idx, key in enumerate(dEa_index):
        rxn = reactionlist[key]
        rxn.CorrectEa = sol_Ea[idx]
    
def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False  
        

if __name__ == '__main__':
    S1 = Species('CO2', 'surface', 'Cu(111)')
    S2 = Species('CO2', 'gaseous')
    S0 = Species('', 'surface')
    R1 = Reaction([(S2,1), (S0,1)], [(S1,1)])
    print(R1)

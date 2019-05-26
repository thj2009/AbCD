import numpy as np

from AbCD.utils import Constant as _const
import AbCD.calculator as _cal

class Reaction(object):
    '''
    Class Reaction
    represent the elementary reaction in surface science
    The reaction kinetics are calculated under the same tempperature,
    Dynamic tempperature Add later
    '''
    def __init__(self):
        # Note: Add name in the input file
        self.name = ''              # name of reaction
        self.reactant = []          # list of reactant
        self.product = []           # list of product
        self.kinetic = {'type': '', 'data': {}}      # Dictionary including kinetic type and parameter
        self.dft_data = {'omega': 0, 'deltaE': 0, 'Ef': 0, 'prefactor': 0}
        self.dE = 0

    def __repr__(self):
        rxn = ''
        for i, react in enumerate(self.reactant):
            stoi = react[1]
            if stoi == 1:
                stoi_spe = str(react[0])
            else:
                stoi_spe = str(stoi) + str(react[0])
            if i == 0:
                rxn += stoi_spe
            else:
                rxn += '+' + stoi_spe
        rxn += '>>'
        for i, prod in enumerate(self.product):
            stoi = prod[1]
            if stoi == 1:
                stoi_spe = str(prod[0])
            else:
                stoi_spe = str(stoi) + str(prod[0])
            if i == 0:
                rxn += stoi_spe
            else:
                rxn += '+' + stoi_spe
        return rxn

    def Arrhenius(self, temp=273.15):
        Arr = {}
        if self.kinetic['type'] == 'Arr':
            Arr['A'] = self.kinetic['data']['A']
            Arr['Ea'] = self.kinetic['data']['Ea']
            Arr['n'] = self.kinetic['data']['n']
        elif self.kinetic['type'] == 'Shomate':
            H_R = self.IS_Enthalpy(temp)
            S_R = self.IS_Entropy(temp)
            H_TS = self.TS_Enthalpy(temp)
            S_TS = self.TS_Entropy(temp)
            Arr['A'] = _const.kb / _const.h * np.exp((S_TS - S_R) / _const.Rg)
            Arr['Ea'] = H_TS - H_R
            Arr['n'] = 1
        return Arr

    def TS_Enthalpy(self, temp=273.15, corr=False):
        if self.kinetic['type'] == 'Arr':
            Hts = self.IS_Enthalpy(temp, corr)
            Hts += self.kinetic['data']['Ea']
            if corr:
                Hts += self.dE
        elif self.kinetic['type'] == 'Shomate':
            if corr:
                Ea0 = _cal.Enthalpy(self.kinetic, temp) - self.IS_Enthalpy(temp)
                Hts = self.IS_Enthalpy(temp, corr)
                Hts += Ea0 + self.dE
            else:
                Hts = _cal.Enthalpy(self.kinetic, temp)
        return Hts

    def TS_Entropy(self, temp=273.15):
        '''
        Calculate the entropy of transition state using shomate parameter
        '''
        Sts = 0
        if bool(self.kinetic['data']):
            Sts = _cal.Entropy(self.kinetic, temp)
        return Sts

    def IS_Enthalpy(self, temp=273.15, corr=False):
        '''
        Calculate the enthalpy of Initial state using shomate parameter
        '''
        His = 0
        for item in self.reactant:
            His += item[1]*item[0].Enthalpy(temp, corr)
        return His

    def IS_Entropy(self, temp=273.15):
        '''
        Calculate the enthalpy of Initial state using shomate parameter
        '''
        Sis = 0
        for item in self.reactant:
            Sis += item[1]*item[0].Entropy(temp)
        return Sis

    def FS_Enthalpy(self, temp=273.15, corr=False):
        '''
        Calculate the entropy of Final state using shomate parameter
        '''
        Hfs = 0
        for item in self.product:
            Hfs += item[1]*item[0].Enthalpy(temp, corr)
        return Hfs

    def FS_Entropy(self, temp=273.15):
        '''
        Calculate the entropy of Final state using shomate parameter
        '''
        Sfs = 0
        for item in self.product:
            Sfs += item[1]*item[0].Entropy(temp)
        return Sfs

    def dH_Enthalpy(self, temp=273.15, corr=False):
        return self.FS_Enthalpy(temp, corr) - self.IS_Enthalpy(temp, corr)

    def dS_Entropy(self, temp=273.15):
        return self.FS_Entropy(temp) - self.IS_Entropy(temp)

    def dG_GibbsFreeEnergy(self, temp=273.15, corr=False):
        dH = self.dH_Enthalpy(temp, corr)
        dS = self.dS_Entropy(temp)
        dG = dH - temp / 1000. * dS
        return dG
    
    def deltaH(self):
        dH = 0
        for item in self.product:
            dH += item[1]*item[0].dE
        for item in self.reactant:
            dH -= item[1]*item[0].dE
        return dH
    
    def deltaEa(self, Tem=298.15, dirc=1):
        if dirc == 1:
            dEa = self.TS_Enthalpy(Tem, True) - self.IS_Enthalpy(Tem, True)
        elif dirc == -1:
            dEa = self.TS_Enthalpy(Tem, True) - self.FS_Enthalpy(Tem, True)
        return dEa
    

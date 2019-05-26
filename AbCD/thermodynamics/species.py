import numpy as np

from AbCD.utils import is_number
from AbCD.utils import Constant as _const
import AbCD.calculator as _cal

class Species(object):
    def __init__(self, name=''):
        self.name = name
        self.formula = ''
        self.thermo = {'type':'', 'data': {}}
        self.vibfreq = []               # vibrational frequency
        self._element = {}
        self.mass = 0                   # mass of speceis unit: g/mol
        self.dE = 0

    def _read_element(self):
        '''
        return the element list given formula
        '''
        spe = str(self.formula)
        element = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0}
        if spe == 'Inert':
            return element
        for idx, ss in enumerate(spe):
            if not is_number(str(ss)):
                stoi = 1
                if not idx == len(spe) - 1:
                    if is_number(spe[idx+1]):
                        stoi = int(spe[idx+1])
                element[ss] += stoi
        self._element = element

    def Enthalpy(self, temp=273.15, corr=False):
        '''
        Calculate Enthapy or heat of formation at specific tempperature
        Unit: kJ/mol
        '''
        H = 0
        if bool(self.thermo['data']):
            H = _cal.Enthalpy(self.thermo, temp)
            if corr:
                H += self.dE
        return H

    def Entropy(self, temp=273.15):
        S = 0
        if bool(self.thermo['data']):
            S = _cal.Entropy(self.thermo, temp)
        return S

    def HeatCapacity(self, temp=273.15):
        Cp = 0
        if bool(self.thermo['data']):
            Cp = _cal.HeatCapacity(self.thermo, temp)
        return Cp

    def ZPEC(self):
        '''
        Calculate Zero point energy correction
        '''
        ZPE = 0
        for v in self.vibfreq:
            nablda = 0.01 / v
            ZPE += 1/2. * _const.h * _const.c / nablda * _const.NA / 1000
        return ZPE

    def Entropy_vib(self, temp=273.15):
        '''
        Calculate Vibrational Entropy from vibrational frequency
        '''
        Svib = 0
        for v in self.vibfreq:
            nablda = 0.01 / v             # wavelength m
            x = _const.h * _const.c / nablda / (_const.kb * temp)
            Svib += x/(np.exp(x)-1) - np.log(1-np.exp(-x))
        Svib *= _const.Rg                     # J mol-1 K-1
        return Svib

class GasSpecies(Species):
    def __init__(self, name=''):
        Species.__init__(self, name)
        self.phase = 'gaseous'

    def __repr__(self):
        return self.formula + '(g)'

    def unicode_repr(self):
        newformula = r''
        for s in self.formula:
            if s in [str(i) for i in range(10)]:
                newformula += '$_' + s + '$'
            else:
                newformula += s
        return newformula + '(g)'

class SurfaceSpecies(Species):
    def __init__(self, name=''):
        Species.__init__(self, name)
        self.phase = 'surface'
        self.site = None
        self.denticity = 0              # Number of site the surface species occupied
        self.bind_info = None            # occupied site geometry configuration and binding energy

    def __repr__(self):
        return self.formula + self.denticity * '*'

    def unicode_repr(self):
        newformula = r''
        for s in self.formula:
            if s in [str(i) for i in range(10)]:
                newformula += '$_' + s + '$'
            else:
                newformula += s
        return newformula + self.denticity * '*'

    def Entropy_2D(self, temp=273.15):
        '''
        Calculate 2D gas Vibrational Entropy
        '''
        S2D = 0
        if self.mass != 0:
            # FIXME: use class site attribute area to calculate
            aa = self.site.lattice_constant['a']
            ll = aa * np.sqrt(2)/2
            Area = np.sqrt(3)/4 * ll **2 * 2
            Strans = _const.Rg * (np.log((2*np.pi*self.mass/1000/_const.NA*_const.kb*temp)/_const.h**2) + np.log(Area) + 2)
            Svib = 0
            for v in self.vibfreq[0:-2]:
                nablda = 0.01 / v      # wavelength m
                x = _const.h * _const.c / nablda / (_const.kb * temp)
                Svib += x/(np.exp(x)-1) - np.log(1-np.exp(-x))
            Svib *= _const.Rg                     # J mol-1 K-1
            S2D = Strans + Svib
        else:
            S2D = 0
        return S2D

    def collisionTheory(self):
        A = 1. / np.sqrt(2 * np.pi * self.mass * _const.u_2_kg * _const.kb) * self.site.area() * 101325
        return A

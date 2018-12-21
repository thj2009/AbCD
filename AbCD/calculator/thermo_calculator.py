import numpy as np
from ..utils import Constant as _const

def Enthalpy(thermo_data, temp=273.15):
    '''
    Calculate enthalpy given temperature and thermodynamic parameter
    '''
    H = 0
    if thermo_data['type'] == 'Shomate':
        tt = temp/1000.
        shomate = thermo_data['data']
        H = shomate['A']*tt + shomate['B']*tt**2/2. + shomate['C']*tt**3/3. + \
            shomate['D']*tt**4/4. - shomate['E']/tt + shomate['F'] - \
            shomate['H'] + shomate['dH']
    elif thermo_data['type'] == 'NASApoly':
        nasapoly = thermo_data['data']
        H = nasapoly['a1'] + nasapoly['a2']*temp/2. + \
            nasapoly['a3']*temp**2/3. + nasapoly['a4']*temp**3/4. + \
            nasapoly['a5']*temp**4/5. + nasapoly['a6']/temp
        H *= _const.Rg * temp/1000
    return H

def Entropy(thermo_data, temp=273.15):
    '''
    Calculate entropy given temperature and thermodynamic parameter
    '''
    S = 0
    if thermo_data['type'] == 'Shomate':
        shomate = thermo_data['data']
        tt = temp/1000.
        S = shomate['A']*np.log(tt) + shomate['B']*tt + shomate['C']*tt**2/2 + \
            shomate['D']*tt**3/3 - shomate['E']/(2*tt**2) + shomate['G']
    elif thermo_data['type'] == 'NASApoly':
        nasapoly = thermo_data['data']
        S = nasapoly['a1'] * np.log(temp) + nasapoly['a2']*temp + \
            nasapoly['a3']*temp**2/2. + nasapoly['a4']*temp**3/3. + \
            nasapoly['a5']*temp**4/4. + nasapoly['a7']
        S *= _const.Rg
    return S

def HeatCapacity(thermo_data, temp=273.15):
    '''
    Calculate heat capacity given temperature and thermodynamic parameter
    '''
    Cp = 0
    if thermo_data['type'] == 'Shomate':
        shomate = thermo_data['data']
        tt = temp/1000.
        Cp = shomate['A'] + shomate['B']*tt + shomate['C']*tt**2 + \
            shomate['D']*tt**3 - shomate['E']/(tt**2)
    elif thermo_data['type'] == 'NASApoly':
        nasapoly = thermo_data['data']
        Cp = nasapoly['a1'] + nasapoly['a2']*temp + nasapoly['a3']*temp**2 \
            + nasapoly['a4']*temp**3 + nasapoly['a5']*temp**4
        Cp *= _const.Rg
    return Cp

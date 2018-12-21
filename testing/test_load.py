import pytest
import numpy as np
import casadi as cas
from AbCD.io_data import In_data
from AbCD.model import CSTR, SimpleKinetic, ReactionNet
from AbCD.visual import parity_plot
import sys
sys.path.append('../')

def test_load():
    sitelist, specieslist, reactionlist = In_data.load_mkm('./WGSdata/')
    assert str(sitelist[0]) == 'Cu-FCC-111-all'
    spenamelist = ['H2(g)', 'CO(g)', 'CO2(g)', 'H2O(g)', 'H*', 'O*', 'OH*', 'H2O*', 'CO*', 'CO2*', 'HCOO*', 'COOH*', '*']
    for spe, spename in zip(specieslist, spenamelist):
        assert str(spe) == spename

    reactionnamelist = [
                'CO(g)+*>>CO*',
                'H2(g)+2*>>2H*',
                'H2O(g)+*>>H2O*',
                'CO2(g)+*>>CO2*',
                'H2O*+*>>H*+OH*',
                'OH*+*>>H*+O*',
                'CO*+O*>>CO2*+*',
                '2OH*>>O*+H2O*',
                'CO2*+H*>>HCOO*+*',
                'CO*+OH*>>COOH*+*',
                'COOH*+*>>CO2*+H*',
                'COOH*+OH*>>CO2*+H2O*'
                ]

    for rxn, rxnname in zip(reactionlist, reactionnamelist):
        assert str(rxn) == rxnname

@pytest.fixture
def cstr():
    sitelist, specieslist, reactionlist = In_data.load_mkm('./WGSdata/')
    conditionlist = In_data.load_condition('./WGSdata/')
    dEa_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    dBE_index = [4, 5 ,6, 7, 8, 9, 10, 11]
    a = CSTR(specieslist=specieslist,
        reactionlist=reactionlist,
        dEa_index=dEa_index,
        dBE_index=dBE_index)
    a.initialize()
    return a

def test_thermoconsis(cstr):
    dE_start = np.zeros(20)
    assert cstr.CheckThermoConsis(dE_start) == False
    dE_corr = cstr.CorrThermoConsis(dE_start)
    assert cstr.CheckThermoConsis(dE_corr) == True

def test_cstr():
    sitelist, specieslist, reactionlist = In_data.load_mkm('./WGSdata/')
    conditionlist = In_data.load_condition('./WGSdata/')
    dEa_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    dBE_index = [4, 5 ,6, 7, 8, 9, 10, 11]
    a = CSTR(specieslist=specieslist,
        reactionlist=reactionlist,
        dEa_index=dEa_index,
        dBE_index=dBE_index)
    a.initialize()
    dE_start = np.zeros(20)
    for condi in conditionlist:
        print(a.fwd_simulation(dE_start, condi))
    expr = []
    mkm = []
    tem = []
    for kk, condi in enumerate(conditionlist):
        mkm.append(a.fwd_simulation(dE_start, condi)['CO(g)'])
        expr.append(condi.TurnOverFrequency['CO(g)'])
#        mkm.append(-condi.SimulatedTOF['H2(g)'])
        tem.append(condi.Temperature)
    
    parity_plot(expr, mkm, 'Experimental TOF', 'Microkinetic TOF')
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(tem, np.log10(mkm), 'k^')
    plt.plot(tem, np.log10(expr), 'ro')

def test_cstr_PE():
    sitelist, specieslist, reactionlist = In_data.load_mkm('./WGSdata/')
    conditionlist = In_data.load_condition('./WGSdata/')
    dEa_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    dBE_index = [4, 5 ,6, 7, 8, 9, 10, 11]
    a = CSTR(specieslist=specieslist,
        reactionlist=reactionlist,
        dEa_index=dEa_index,
        dBE_index=dBE_index)
    a.initialize()
    dE_start = np.random.uniform(-20, 20, 20)
    dE_start = a.CorrThermoConsis(dE_start)
#    print(dE_start)
#    print(a.CheckThermoConsis(dE_start))
#    for condi in conditionlist:
#        print(a.parameter_estimation(dE_start, conditionlist, print_screen=True))
    dE, obj = a.parameter_estimation(dE_start, conditionlist, print_screen=True)
    print(dE, obj)
    # parity Plot
    expr = []
    mkm = []
    for kk, condi in enumerate(conditionlist):
        mkm.append(a.fwd_simulation(dE, condi)['CO(g)'])
        expr.append(condi.TurnOverFrequency['CO(g)'])
#        mkm.append(-condi.SimulatedTOF['H2(g)'])
    
    parity_plot(expr, mkm, 'Experimental TOF', 'Microkinetic TOF')

def test_cstr_MS_PE():
    sitelist, specieslist, reactionlist = In_data.load_mkm('./WGSdata/')
    conditionlist = In_data.load_condition('./WGSdata/')
    dEa_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    dBE_index = [4, 5 ,6, 7, 8, 9, 10, 11]
    a = CSTR(specieslist=specieslist,
        reactionlist=reactionlist,
        dEa_index=dEa_index,
        dBE_index=dBE_index)
    a.initialize(scale=1e-4)
#    dE_start = np.zeros(20)
#    a.parameter_estimation(dE_start, conditionlist, print_screen=True, report='%d.txt'%0)

    for i in range(200):
        dE_start = np.random.uniform(-20, 20, 20)
        dE_start = a.CorrThermoConsis(dE_start)
        a.parameter_estimation(dE_start, conditionlist, print_screen=True, report='./multi-start/%d.txt'%i)


def test_new_pe():
    sitelist, specieslist, reactionlist = In_data.load_mkm('./WGSdata/')
    conditionlist = In_data.load_condition('./WGSdata/')
    dEa_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    dBE_index = [4, 5 ,6, 7, 8, 9, 10, 11]
    a = CSTR(specieslist=specieslist,
        reactionlist=reactionlist,
        dEa_index=dEa_index,
        dBE_index=dBE_index)
    a.initialize()
    dE_start = np.random.uniform(-20, 20, 20)
    dE_start = a.CorrThermoConsis(dE_start)

    evidence = {}
    evidence['type'] = 'rel'
    evidence['err'] = 1
    
    prior = {}
    prior['type'] = 'Ridge'
    prior['L2'] = 1e-5
    prior['lbound'] = [-20] * 20
    prior['ubound'] = [20] * 20

    dE, obj = a.mle_estimation(dE_start, conditionlist, evidence, prior)
    expr = []
    mkm = []
    for kk, condi in enumerate(conditionlist):
        mkm.append(a.fwd_simulation(dE, condi)['CO(g)'])
        expr.append(condi.TurnOverFrequency['CO(g)'])
#        mkm.append(-condi.SimulatedTOF['H2(g)'])
    
    parity_plot(expr, mkm, 'Experimental TOF', 'Microkinetic TOF')


if __name__ == '__main__':
#    test_cstr_MS_PE()
#    test_cstr()
#    test_cstr_PE()
    test_new_pe()

# test reaction reaction network
import pytest
import sys
sys.path.append('../')

import AbCD
from AbCD import SimpleKinetic, CSTR
import numpy as np

@pytest.fixture
def network():

    CO = AbCD.SurfaceSpecies()
    CO.formula = 'CO'
    CO.denticity = 1
    CO.thermo = {'type':'Shomate', 'data':
        {'dH':-57015.02, 'dS': 202.05, 'A': 27.01916, 'B': 0.39882, 'C': 10.80028,
         'D':-5.39747, 'E': 0.11318, 'F': -57077.85, 'G':234.83, 'H':-57015.02}}
    O = AbCD.SurfaceSpecies()
    O.formula = 'O'
    O.denticity = 1
    O.thermo = {'type':'Shomate', 'data':
        { 'dH':-42119.61, 'dS':13.85, 'A':18.66427, 'B': 18.08991, 'C': -19.15144,
         'D': 7.00167, 'E': -0.29627, 'F': -42128.91229, 'G': 30.16, 'H': -42121.71}}

    CO2 = AbCD.SurfaceSpecies()
    CO2.formula = 'CO2'
    CO2.denticity = 1
    CO2.thermo = {'type':'Shomate', 'data':
        {'dH':-99310.40, 'dS':218.57, 'A':27.72336, 'B': 48.60556, 'C':-28.93575,
         'D':6.74397, 'E':-0.18280, 'F':-99321.19, 'G':237.82, 'H':-99310.40}}
    R1 = AbCD.Reaction()
    R1.reactant = [(CO, 1), (O, 1)]
    R1.product = [(CO2, 1)]
    R1.kinetic = {'type':'Arr', 'data': {'A': 1e13, 'Ea': 10, 'n': 1.0}}
    specieslist=[CO, O, CO2]
    reactionlist = [R1]
    rxn_net = CSTR(specieslist, reactionlist, [0], [0,1,2])
    return rxn_net

def test_reaction_species_count(network):
    assert network.ngas == 0
    assert network.nsurf == 3
    assert network.nrxn == 1

def test_reaction_parameter_count(network):
    assert network.dEa_index == [0]
    assert network.dBE_index == [0, 1, 2]

def test_stoimat(network):
    network._stoimat_generate()
    np.testing.assert_array_equal(network.stoimat, np.array([[-1, -1, 1]]))

def test_energy_expression(network):
    network.build_kinetic()
    print(network._kf)
    print(network._kr)
    print(network._Keq)

def test_rate_expression(network):
    network.build_kinetic()
    network.build_rate()
    print(network._rate)
    print(network._rfor)

#def test
if __name__ == '__main__':
    n = network()
#    for spe in n.specieslist:
#        print(spe.Enthalpy)
#    test_energy_expression(n)
#    test_rate_expression(n)
    n.initialize()
    print(n._dae_)


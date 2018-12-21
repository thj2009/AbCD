# test case try pytest fixture
import sys
sys.path.append('../')
import pytest
from AbCD import GasSpecies, SurfaceSpecies, Reaction

@pytest.fixture
def empty_Gas():
    return GasSpecies()

@pytest.fixture
def empty_Surf():
    return SurfaceSpecies()

@pytest.fixture
def empty_Rxn():
    return Reaction()

def test_empty_species_output(empty_Gas, empty_Surf):
    assert str(empty_Gas) == '(g)'
    assert str(empty_Surf) == ''

@pytest.mark.parametrize("formula,denticity,expected_surf", [
    ('CO', 0, 'CO'),
    ('O', 1, 'O*'),
])

def test_nonempty_surface_output(empty_Surf, formula, denticity, expected_surf):
    empty_Surf.formula = formula
    empty_Surf.denticity = denticity
    assert str(empty_Surf) == expected_surf

def test_species_energy(empty_Gas):
    assert empty_Gas.Enthalpy() == 0

def test_empty_reaction_output(empty_Rxn):
    assert str(empty_Rxn) == '>>'

def test_nonempty_reaction_output(empty_Rxn):
    R1 = GasSpecies()
    R1.formula = 'CO'

    R2 = SurfaceSpecies()
    R2.formula = 'O'
    R2.denticity = 2

    R3 = GasSpecies()
    R3.formula = 'CO2'

    empty_Rxn.reactant = [(R1, 1), (R2, 1)]
    empty_Rxn.product = [(R3, 1)]
    assert str(empty_Rxn) == 'CO(g)+O**>>CO2(g)'


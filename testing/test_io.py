import pytest 
import tempfile
from AbCD.io_data import In_data
import json


@pytest.fixture
def site():
    site_info = {}
    site_info['name'] = 'Cu'
    site_info['metal'] = 'Cu'
    site_info['struct'] = 'FCC'
    site_info['facet'] = '111'
    site_info['lattice_constant'] = {}
    site_info['mass'] = 64
    site_string = json.dumps([site_info])
    return site_string
   
@pytest.fixture
def species():
    species_info = {}
    species_info['name'] = 'carbonmonoxide'
    species_info['formula'] = 'CO'
    species_info['phase'] = 'g'
    species_info['mass'] = 0
    species_info['site_name'] = 'Cu'
    species_string = json.dumps([species_info])
    return species_string

def test_read_site(site):
    sitelist = In_data.read_site(site)
    print(str(sitelist[0]))
    assert len(sitelist) == 1
    assert sitelist[0].name == 'Cu'
    assert sitelist[0].metal == 'Cu'
    assert sitelist[0].struct == 'FCC'
    assert sitelist[0].facet == '111'
    
def test_read_species(species, site):
    sitelist = In_data.read_site(site)
    specieslist = In_data.read_species(species, sitelist)
    assert len(specieslist) == 2
    assert specieslist[0].name == 'carbonmonoxide'
    assert str(specieslist[0]) == 'CO(g)'
    assert specieslist[0].phase == 'gaseous'
    
    assert specieslist[1].name == '*'
    assert str(specieslist[1]) == '*'
    assert specieslist[1].phase == 'surface'
    


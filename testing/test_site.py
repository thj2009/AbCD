# test active site and species
import sys
sys.path.append('../')
import AbCD

Cu = AbCD.ActiveSite()

def test_activesite():
    Cu.name = 'Cu'
    Cu.struct = 'FCC'
    Cu.metal = 'Cu'
    Cu.facet = '111'
    Cu.site = '3fold'
    assert str(Cu) == 'Cu-FCC-111-3fold'




def test_lattice():
    top = AbCD.ActiveSite('top')
    Latt = AbCD.Lattice()
    Latt.cell_vector = [[1, 0], [0, 1]]
    Latt.repeat = [10, 10]
    
    Latt.sites = [top]
    Latt.unit_cell = {'top': [0, 0]}
    Latt.neighbs = [('top/00', 'top/01'), ('top/00', 'top/10)')]
    
    Latt.initialization()
#    print(Latt.occupancy)
if __name__ == '__main__':
    test_lattice()
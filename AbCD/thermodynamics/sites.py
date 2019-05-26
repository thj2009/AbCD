import numpy as np

class ActiveSite(object):
    '''
    Active Site information
    '''
    def __init__(self, name=''):
        # HINT: may need to generalize to other typer of active site
        # HINT: create inherenty class from this site class
        self.name = name
        self.metal = ''
        self.struct = ''
        self.facet = ''
        self.site = ''
        self.lattice_constant = {}
        self.mass = 0

    def __str__(self):
        return '-'.join([self.metal, self.struct, self.facet, self.site])

    def detail_repr(self):
        '''
        detail active site representation
        '''
        # TODO add this method
        pass

    def area(self):
        '''
        Area per site, used for 2D translation entropy calculation
        '''
        a = self.lattice_constant['a']
        d = a / np.sqrt(2)
        area = np.sqrt(3) / 2 * d**2
        return area

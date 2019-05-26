from AbCD.io_data import In_data
from AbCD.utils import Constant as _c
import numpy as np

sitelist, specieslist, reactionlist = In_data.load_mkm('./KineticExp/MetSyndata/')
conditionlist = In_data.load_condition('./KineticExp/MetSyndata/')

aa = 18785400000
dev = 9.87845476211e+13
for spe in specieslist:
    if spe.phase == 'gaseous':
        A = 1 / np.sqrt(2 * np.pi * spe.mass /1000./_c.NA * _c.kb)
        print(spe, A/dev)
#        if spe.name == 'H2(g)':
#            print(A)
#            print(A/aa)


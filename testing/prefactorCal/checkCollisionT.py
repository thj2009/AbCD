# %%
#KINETIC MODEL


from AbCD.model import CSTR
from AbCD.io_data import In_data
from AbCD.visual import tpd_profile
import numpy as np

from AbCD.model import CSTR
from AbCD.io_data import In_data
from AbCD.visual import tpd_profile
import numpy as np

sitelist, specieslist, reactionlist = In_data.load_mkm('./KineticExp/C1Kinetic/')

for spe in specieslist:
    if spe.phase != 'gaseous':
        A = spe.collisionTheory()
        print(str(spe) + '%.6e' %A)
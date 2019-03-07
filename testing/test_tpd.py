from AbCD.model import VacuumTPD
from AbCD.io_data import In_data

sitelist, specieslist, reactionlist = In_data.load_mkm('./TPDdata/')
conditionlist = In_data.load_condition('./TPDdata/')

dEa_index = []
dBE_index = [0]
a = VacuumTPD(specieslist=specieslist,
              reactionlist=reactionlist,
              dEa_index=dEa_index,
              dBE_index=dBE_index)
a.initialize()

dP = [0]
a.fwd_simulation(dP, conditionlist[0])
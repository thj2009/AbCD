from AbCD.model import VacuumTPD
from AbCD.io_data import In_data
from AbCD.visual import tpd_profile
sitelist, specieslist, reactionlist = In_data.load_mkm('./testing/TPDdata/')
conditionlist = In_data.load_condition('./testing/TPDdata/')

dEa_index = []
dBE_index = [1]
a = VacuumTPD(specieslist=specieslist,
              reactionlist=reactionlist,
              dEa_index=dEa_index,
              dBE_index=dBE_index)
a.initialize(pump_level=1)

dP = [-10]
# out = a.fwd_simulation(dP, conditionlist[0])

# _, nd = out.shape
# des = [2 * (-out[1, i+1] + out[1, i]) for i in range(nd - 1)]
#
# tpd_profile(conditionlist[0].TemGrid[:-1], des, r_=[100, 500])
import pytest
from AbCD.io_data import In_data
from AbCD import CSTR
import numpy as np

sitelist, specieslist, reactionlist = In_data.load_mkm('./testing/WGSdata/')
conditionlist = In_data.load_condition('./testing/WGSdata/')

dBE_index = [4, 5, 6, 7, 8, 9, 10, 11]
dEa_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]


def test_cstr_build():
    wgs_cstr = CSTR(
        specieslist=specieslist,
        reactionlist=reactionlist,
        dBE_index=dBE_index,
        dEa_index=dEa_index
    )
    assert wgs_cstr.Pnlp.shape[0] == len(dEa_index) + len(dBE_index)
    wgs_cstr.initialize()

    E = 0 * np.ones(len(dEa_index) + len(dBE_index))
    E = wgs_cstr.CorrThermoConsis(E)
    assert wgs_cstr.CheckThermoConsis(E)

    for condi in conditionlist:
        print wgs_cstr.fwd_simulation(E, condi, reltol=1e-12, abstol=1e-14)






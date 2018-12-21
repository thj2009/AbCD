# -*- coding: utf-8 -*-
"""
Microkinetic Modeling 
Water gas shift application
Input:
    Reaction network
    DFT data
    Reactor data
Output:
    General MKM model,
    thermal dynamic consistency constraint
    Link to all thermal data
@author: hut216
"""

import sys
sys.path.append(r'D:\Anaconda2\Lib\site-packages\casadi_2.4.2')
#sys.path.append('/home/hut216/CASADI/casadi_2.4.2')

import casadi as cas
import MF_MKM
import read_file
import numpy as np

# Physical Constant
R_g = 8.314
NA = 6.02e23

# Define
speciesfile = './WGSdata/WGS_species.txt'
reactionfile = './WGSdata/WGS_rxn_TS.txt'

Spe_list = read_file.import_species(speciesfile)
Rxn_list = read_file.import_reaction(reactionfile, Spe_list)

WGS = MF_MKM.ReactionNet(Spe_list, Rxn_list)
WGS.StoiMat()
#WGS.printscreen()


# Define Experiment Parameter
inlet_flow = 118./1e6/60.  # m^3/s 
inlet_flow *= 10000   # Change the flow rate, make sure the rate is initial rate
Tem = 190 + 273.15  # Kelvin
inlet_press = 1  # atm
inlet_mole = inlet_flow * inlet_press * 1.01325e5/(R_g*Tem)       # mole/s

cata_mass = 0.2  # g
cata_per_area = 9.6  # m^2/g
cata_density = 1.8e19  # atoms/m^2
cata_mole = cata_mass*cata_per_area*cata_density/NA  # mole

reactor = dict()
reactor['inletflow'] = inlet_flow
reactor['inletmole'] = inlet_mole
reactor['inletpress'] = inlet_press
reactor['temperature'] = Tem
reactor['cata_mole'] = cata_mole

 
dEa_index = [4, 5, 6, 7, 8, 9, 10, 11]
dBE_index = [4, 5, 6, 7, 8, 9, 10, 11]

WGS_CSTR = MF_MKM.CSTR_model(WGS, reactor, dEa_index, dBE_index)
WGS_CSTR.Build_model()

experimetnfile = './WGSdata/Ribeiro_EXP_data.txt'
Ppress_in_m, Tor_m = read_file.import_experiment(experimetnfile)
Ppress_in_m = np.array(Ppress_in_m).T

#AA = [2.7192734440845436, 0.6208656139022156, 4.773293537215737, 8.478309093225116, 1.046958753093239, -3.509453123089648, 5.3783819859486846, 10.811175528055154, -3.439886563182078, -7.055636227902168, 4.925109892574479, -18.250643941541668, -5.212998536415709, -6.510623434165154, -2.0897520407017223, 1.4206542742160804];
AA = np.zeros(16).tolist()
AA = [-6.41689051525443, 8.056892120140944, -8.043707817862952, -6.446591396561113, -10.681555250261903, -18.583086296251174, -13.22788041045843, 7.520732322254317, 7.041109331765269, -10.570224639411466, 1.6392485238501802, -9.598115960715935, -17.335879158811867, 10.700342334683317, -19.412096158618727, -15.452219766792501]
dEa_ = AA[0:8]
dBE_ = AA[8:]
Tor_s = []
for k in range(19):
    P_flow = Ppress_in_m[:,k]; 
    print(P_flow)
#    dEa_, dBE_ = WGS_CSTR.CorrThermoConsis(dEa_, dBE_, fix_BE=dBE_index)
    print(WGS_CSTR.ForwardCal(dEa_, dBE_, P_flow))
    Tor_s.append(WGS_CSTR.ForwardCal(dEa_, dBE_, P_flow))
    
    
#Tor_s0 = [0.2125406078787925, 0.15562335295230292, 0.24011235785655402, 0.19039216918932006, 0.26530192735072677, 0.19447804297105623, 0.19994595639523155, 0.19738947906547702, 0.19987469968532875, 0.19963840210935319, 0.21842108981534758, 0.2076799693901096, 0.18479378970858584, 0.20237739195653665, 0.076175481765535957, 0.32477633462863836, 0.11496759844012927, 0.23461007572021605, 0.15820816498060936];
Tor_s0 = [10.560558907432801, 10.227254522272148, 10.607219113897219, 10.449124236142289, 10.572338659595582, 10.474199692087483, 10.281549229322616, 10.656832021880495, 10.42507292925162, 10.523448479470447, 10.913735286400566, 10.769235314092651, 10.279922841250688, 10.635655803289982, 3.8188623252118341, 18.795554269291422, 5.8801033825679694, 12.918338854380393, 8.3426550251757359];
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
Tor_m = np.array(Tor_m)
#plt.figure()
fig = plt.figure()
ax = plt.gca()
plt.plot(Tor_m, Tor_s, 'bo')
d = Tor_m.max() - Tor_m.min()
plt.plot([Tor_m.min() - 0.1*d, Tor_m.max() + 0.1*d], [Tor_m.min() - 0.1*d, Tor_m.max() + 0.1*d], 'k:', linewidth = 2)
plt.xlim([Tor_m.min() - 0.1*d, Tor_m.max() + 0.1*d])
#    plt.title('Median TOF', fontsize = 14)
plt.xlabel(r'Experimental TOF (10$^{-2}$ s$^{-1}$)', fontsize=16)
plt.ylabel(r'Calculated TOF (10$^{-2}$ s$^{-1}$)', fontsize=16)
plt.legend(fontsize=12)
# this is an inset axes over the main axes
inset_axes = inset_axes(ax, 
                    width="50%", # width = 30% of parent_bbox
                    height=2.0, # height : 1 inch
                    loc=2)
plt.plot(Tor_m, Tor_s0, 'ro')
plt.plot([min(Tor_m)-0.1*d, max(Tor_m)+0.1*d], [min(Tor_m)-0.1*d, max(Tor_m)+0.1*d], 'k:')
plt.xlim([min(Tor_m)-0.1*d, max(Tor_m)+0.1*d])
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.xticks([])
plt.yticks([])

print np.sum((np.array(Tor_m) - np.array(Tor_s))**2)
print np.exp(abs(np.log(Tor_m) - np.log(Tor_s)).mean())-1




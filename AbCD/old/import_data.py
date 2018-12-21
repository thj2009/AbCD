# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:31:03 2017

@author: huijie
"""

# import species data
import numpy as np

def import_data():
    
    spe_list =[]
    name_list = []
#    spe_list.append(('Arg', {})) # Add inertial gas in the species list
#    name_list.append('Arg')
    
    f1 = open(".\data\WGS_species.txt", 'r')
#    f1 = open("./data/WGS_species.txt", 'r')
    s1 = f1.readlines()
    for v1 in s1:
        sh = v1.split()
        shomate = dict()
        if sh !=[]:
            shomate['A'] = float(sh[sh.index('A')+1])
            shomate['B'] = float(sh[sh.index('B')+1])
            shomate['C'] = float(sh[sh.index('C')+1])
            shomate['D'] = float(sh[sh.index('D')+1])
            shomate['E'] = float(sh[sh.index('E')+1])
            shomate['F'] = float(sh[sh.index('F')+1])
            shomate['G'] = float(sh[sh.index('G')+1])
            shomate['H'] = float(sh[sh.index('H')+1])
            shomate['dH'] = float(sh[sh.index('dH')+1])
            S = (sh[0], shomate)
            spe_list.append(S)
            name_list.append(sh[0])
    f1.close()
    
    spe_list.append(('*',{}))   # Add clean surface in the species list
    name_list.append('*')
    
    # import reaction data and return kinetic data

    f2 = open(".\data\WGS_rxn.txt", 'r')
#    f2 = open("./data/WGS_rxn.txt", 'r')
    s2 = f2.readlines()
    kine_list = []
    for v2 in s2:
        ss = v2.split()
        ii = ss.index('A')
        kin = dict()
        kin['A'] = float(ss[ss.index('A')+1])
        kin['Ea'] = float(ss[ss.index('Ea')+1])
        kin['n'] = float(ss[ss.index('n')+1])
        kine_list.append(kin)
    
    f3 = open('.\data\stoichiometric.txt', 'r')
#    f3 = open('./data/stoichiometric.txt', 'r')
    s3 = f3.readlines()
    stoi = np.zeros((len(s2), len(spe_list)))
    for i in range(len(s3)):
        ss = s3[i].split()
#        print(ss)
        sss = [int(ii) for ii in ss]
#        print(sss)
        stoi[i] = sss
        
#    f4 = open('.\data\Ribeiro_EXP_data_Arbalance.txt', 'r')
    f4 = open('.\data\Ribeiro_EXP_data.txt', 'r')
#    f4 = open('./data/Ribeiro_EXP_data.txt', 'r')
    s4 = f4.readlines()
    k_exp = len(s4)
    j_g = 4
    flowin_m = []
    tor = []
    for i in range(k_exp):
        ss = s4[i].split()
        sss = [float(ii) for ii in ss]
        flowin_m.append(sss[0:-1])
        tor.append(sss[-1])

    return(spe_list, name_list,kine_list,stoi, flowin_m, tor)

    
if __name__ == '__main__':
    (a,b,c,d,e,f) = import_data()
    print(a)
    print(b)
    print(c)
    print(d)
    print(e)
    print(f)
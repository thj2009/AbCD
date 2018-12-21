# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:31:03 2017

@author: huijie
"""

# import species data
import numpy as np
import MF_MKM


def is_number(s):
#    print(s)
    try:
        int(s)
        return True
    except ValueError:
#        raise
        return False
        
def import_species(speciesfile):
    
    spe_list =[]
    
    file = open(speciesfile, 'r')
    lines = file.readlines()
    for line in lines:
        sline = line.split()
        shomate = dict()
        if sline !=[]:
            
            shomate['A'] = float(sline[sline.index('A')+1])
            shomate['B'] = float(sline[sline.index('B')+1])
            shomate['C'] = float(sline[sline.index('C')+1])
            shomate['D'] = float(sline[sline.index('D')+1])
            shomate['E'] = float(sline[sline.index('E')+1])
            shomate['F'] = float(sline[sline.index('F')+1])
            shomate['G'] = float(sline[sline.index('G')+1])
            shomate['H'] = float(sline[sline.index('H')+1])
            shomate['dH'] = float(sline[sline.index('dH')+1])
            
            if sline[0][-1] == '*':
                name = sline[0][:-1]
                Spe = MF_MKM.Species(name, 'surf', shomate = shomate)
            elif sline[0][-1] == 'g':
                name = sline[0][:-1]
                Spe = MF_MKM.Species(name, 'g', shomate = shomate)
            spe_list.append(Spe) 
    Spe = MF_MKM.Species('', 'surf')
    spe_list.append(Spe)  
    file.close()
    return spe_list

def import_reaction(reactionfile, specieslist):
    # import reaction data and return kinetic data

    file = open(reactionfile, 'r')
    lines = file.readlines()
    Rxn_list = []
    for line in lines:
        sline = line.split('//')
        if sline !=[] and sline != ['\n']:
#            print sline
            rxnname = sline[0].strip()
            kinline = sline[2].split()
            kinetictype = kinline[0]
            if kinetictype == 'Arr':
                param = dict()
                param['A'] = float(kinline[kinline.index('A')+1])
                param['Ea'] = float(kinline[kinline.index('Ea')+1])
                param['n'] = float(kinline[kinline.index('n')+1])
            elif kinetictype == 'Shomate':
                param = dict()
                param['A'] = float(kinline[kinline.index('A')+1])
                param['B'] = float(kinline[kinline.index('B')+1])
                param['C'] = float(kinline[kinline.index('C')+1])
                param['D'] = float(kinline[kinline.index('D')+1])
                param['E'] = float(kinline[kinline.index('E')+1])
                param['F'] = float(kinline[kinline.index('F')+1])
                param['G'] = float(kinline[kinline.index('G')+1])
                param['H'] = float(kinline[kinline.index('H')+1])
                param['dH'] = float(kinline[kinline.index('dH')+1])
            kin = {}
            kin['type'] = kinetictype
            kin['param'] = param
            rxn = sline[1].split()
            react = rxn[:rxn.index('>>')]; react = ''.join(react).split('+')
            prod = rxn[rxn.index('>>')+1:]; prod = ''.join(prod).split('+')

            Reactant = []
            for rr in react:
                if is_number(rr[0]):
                    spe = rr[1:]; stoi = int(rr[0])
                    index = MF_MKM.get_index_species(spe, specieslist)
                else:
                    spe = rr; stoi = 1
                    index = MF_MKM.get_index_species(spe, specieslist)
                Reactant.append((specieslist[index], stoi))
            Prod = []
            for pp in prod:
                if is_number(pp[0]):
                    spe = pp[1:]; stoi = int(pp[0])
                    index = MF_MKM.get_index_species(spe, specieslist)
                else:
                    spe = pp; stoi = 1
                    index = MF_MKM.get_index_species(spe, specieslist)
                Prod.append((specieslist[index], stoi))
            Rxn = MF_MKM.Reaction(Reactant, Prod, kin)
            Rxn.name = rxnname
            Rxn_list.append(Rxn)
    file.close()
    return Rxn_list

def import_vib(vibfile , specieslist):
    file = open(vibfile, 'r')
    lines = file.readlines()
    for line in lines:
        sline = line.split()
        if sline != []:
            spe = sline[0]
            vib = [float(v) for v in sline[sline.index('//')+1:]]
            index = MF_MKM.get_index_species(spe, specieslist)
            if index != None:
                specieslist[index].vibfreq = vib

def import_rxn_DFT(RxnDFTfile, rxnlist):
    file = open(RxnDFTfile, 'r')
    lines = file.readlines()
    for line in lines:
        sline = line.split('//')
        if sline !=[]:
            rxnname = sline[0].strip()
            dataline = sline[2].split()
            deltaE = float(dataline[dataline.index('deltaE')+1]) * 96.486
            Ef = float(dataline[dataline.index('Ef')+1]) * 96.486
            omega = float(dataline[dataline.index('omega')+1])
            AA = float(dataline[dataline.index('A')+1])
            idx = MF_MKM.get_index_reactions(rxnname, rxnlist)
            rxnlist[idx].deltaE = deltaE
            rxnlist[idx].Ef = Ef
            rxnlist[idx].omega = omega
            rxnlist[idx].prefactor = AA
                
def import_BE(BEfile , specieslist):
    file = open(BEfile, 'r')
    lines = file.readlines()
    for line in lines:
        sline = line.split()
        if sline != []:
            spe = sline[0]
            if sline[sline.index('//')+1] == 'BEev':
                index = MF_MKM.get_index_species(spe, specieslist)
                specieslist[index].BindEnergy = float(sline[sline.index('BEev')+1]) * 96.485
        
def import_experiment(experimentfile):   
    file = open(experimentfile, 'r')
    lines = file.readlines()
    k_exp = len(lines)
    Ppress_in_m = []
    Tor_m = []
    for line in lines:
        sline = line.split()
        sss = [float(ii) for ii in sline]
        Ppress_in_m.append(sss[0:-1])
        Tor_m.append(sss[-1])

    return(Ppress_in_m, Tor_m)
    
def check_comment(line):
    if line[0] == '#':
        return True
    

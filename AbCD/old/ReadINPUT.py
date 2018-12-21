# -*- coding: utf-8 -*-
"""
Created on Tue May 01 21:59:05 2018

read Composited input file 
comment line starting with '#'

section seperation
start 'sectionname'
end  'sectionname'

List of section name
1. activesite
2. species
3. reaction 
4. reactor
5. condition
    5.a inlet
    5.b outlet
6. vibrational frequency of species
7. DFT of reaction 
8. species thermodynamic parameter
9. reaction kinetic parameter
10. bind information
11. parameter
    11.a species parameter
    11.b species parameter
@author: huijie
"""

import MF_MKM
import numpy as np

s_comment = '#'
s_start = 'start'
s_end = 'end'
s_seperate = '//'

def read_file(composefile):
    ActiveSite = MF_MKM.ActiveSite()
    SpeciesList = []
    ReactionList = []
    ConditionList = []
    SpeciesParameter = []
    ReactionParameter = []
    with open(composefile, 'r') as fp:
        flag_activesite = 0
        flag_species = 0        
        flag_reaction = 0
        flag_reactor = 0
        flag_vibration = 0
        flag_rxnDFT = 0
        flag_speciesthermo = 0
        flag_reactionthermo = 0
        flag_condition = 0
        flag_inlet = 0
        flag_outlet = 0
        flag_bindinfo = 0
        flag_parameter = 0
        flag_speciesparameter = 0
        flag_reactionparameter = 0
        count = 0
        for line in fp:
            count += 1
            line = line.strip()
            
            # Check Start Section
            if line.startswith(s_comment):
                continue
            
            # Check End Section
            if line.startswith(s_end):
                # raise exception if start and end not match
                sl = line.split()
                end_section = sl[1]
                if end_section == 'activesite' and flag_activesite == 1:
                    flag_activesite = 0
                elif end_section == 'species' and flag_species == 1:
                    spe = MF_MKM.Species('*', 'surface')
                    SpeciesList.append(spe)
                    flag_species = 0
                elif end_section == 'reaction' and flag_reaction == 1:
                    flag_reaction = 0
                elif end_section == 'reactor' and flag_reactor == 1:
                    flag_reactor = 0
                elif end_section == 'condition' and flag_condition == 1:
                    if len(ConditionList) != num_condition:
                        raise NameError('Number of reaction Condition not Match')
                    flag_condition = 0
                elif end_section == 'vibration' and flag_vibration == 1:
                    flag_vibration = 0
                elif end_section == 'reaction_DFT' and flag_rxnDFT == 1:
                    flag_rxnDFT = 0
                elif end_section == 'reaction_thermo' and flag_reactionthermo == 1:
                    flag_reactionthermo = 0
                elif end_section == 'species_thermo' and flag_speciesthermo == 1:
                    flag_speciesthermo = 0
                elif end_section == 'inlet' and flag_inlet == 1:
                    flag_inlet = 0
                elif end_section == 'outlet' and flag_outlet == 1:
                    flag_outlet = 0
                elif end_section == 'bindinfo' and flag_bindinfo == 1:
                    flag_bindinfo = 0
                elif end_section == 'parameter' and flag_parameter == 1:
                    flag_parameter = 0
                elif end_section == 'species_parameter' and flag_speciesparameter == 1:
                    flag_speciesparameter = 0
                elif end_section == 'reaction_parameter' and flag_reactionparameter == 1:
                    flag_reactionparameter = 0
                else:
                    # Raise error
                    print 'Line %d: Section not match' %count  
                    
                    
            if flag_activesite == 1:
            # extract active site information
                if line.split() == []:
                    continue
                sl = line.split(s_seperate)
                # name of metal // mass // struct // facet // Lattice constant
                ActiveSite.metal = sl[0].replace(' ', '')
                ActiveSite.mass = float(sl[1].replace(' ', ''))
                ActiveSite.struct = sl[2].replace(' ', '')
                ActiveSite.facet = sl[3].replace(' ', '')
                lc = {}
                key = sl[4].split()[0]
#                    print key
                value = float(sl[4].split()[1])
                lc[key] = value
                ActiveSite.LatticeConstant = lc
                
            if flag_species == 1:
            # extract species information
                if line.split() == []:
                    continue
                species = MF_MKM.Species()
                sl = line.split(s_seperate)
                # name of species // mass // phase // denticity
                species.name = sl[0].replace(' ', '')
                species.mass = float(sl[1].replace(' ', ''))
                species.phase = sl[2].replace(' ', '')
                species.denticity = int(sl[3].replace(' ', ''))
                SpeciesList.append(species)
                if species.phase == 'surface':
                    species.site = ActiveSite
#                    print species.site
                species.element = species.readelement()
            
            if flag_reaction == 1:
            # extract reaction information
                if line.split() == []:
                    continue
                reaction = MF_MKM.Reaction()
                sl = line.split(s_seperate)
                # name of reaction // reaction
                # reaction: reactant >> product
                reaction.name = sl[0].replace(' ', '')
                rxn = sl[1].split()
                react = rxn[:rxn.index('>>')]; react = ''.join(react).split('+')
                prod = rxn[rxn.index('>>')+1:]; prod = ''.join(prod).split('+')
                React = []
                for rr in react:
                    if is_number(rr[0]):
                        spe = rr[1:]; stoi = int(rr[0])
                        index = MF_MKM.get_index_species(spe, SpeciesList)
                    else:
                        spe = rr; stoi = 1
                        index = MF_MKM.get_index_species(spe, SpeciesList)
                    React.append((SpeciesList[index], stoi))
                Prod = []
                for pp in prod:
                    if is_number(pp[0]):
                        spe = pp[1:]; stoi = int(pp[0])
                        index = MF_MKM.get_index_species(spe, SpeciesList)
                    else:
                        spe = pp; stoi = 1
                        index = MF_MKM.get_index_species(spe, SpeciesList)
                    Prod.append((SpeciesList[index], stoi))
                    
                reaction.reactant = React
                reaction.product = Prod
                ReactionList.append(reaction)
                  
            if flag_vibration == 1:
            # extract vibration information
                if line.split() == []:
                    continue
                sl = line.split(s_seperate)
                spe = sl[0].replace(' ', '')
                freq = sl[1].split()
                idx = MF_MKM.get_index_species(spe, SpeciesList)
                vib = [float(v) for v in freq]
                SpeciesList[idx].vibfreq = vib
            
            if flag_rxnDFT == 1:
            # extract reaction DFT calculation data
                if line.split() == []:
                    continue
                sl = line.split(s_seperate)
                # name //  deltaE // Ef // omega // prefactor
                rxn = sl[0].replace(' ', '')
                idx = MF_MKM.get_index_reactions(rxn, ReactionList)
                ReactionList[idx].deltaE = float(sl[1].replace(' ', ''))
                ReactionList[idx].Ef = float(sl[2].replace(' ', ''))
                ReactionList[idx].omega = float(sl[3].replace(' ', ''))
                ReactionList[idx].prefactor = float(sl[4].replace(' ', ''))
            
            if flag_speciesthermo == 1:
            # extract species thermodynamic data
                if line.split() == []:
                    continue
                sl = line.split(s_seperate)
                # name of species // type of data // thermodynamic data
                spe = sl[0].replace(' ', '')
                idx = MF_MKM.get_index_species(spe, SpeciesList)
                thermotype = sl[1].replace(' ', '')
                
                thermo = {}
                thermo['type'] = thermotype
                thermo['data'] = {}
                if thermotype == 'Shomate':
                    shomate = {}
                    sline = sl[2].split()
                    shomate['A'] = float(sline[sline.index('A')+1])
                    shomate['B'] = float(sline[sline.index('B')+1])
                    shomate['C'] = float(sline[sline.index('C')+1])
                    shomate['D'] = float(sline[sline.index('D')+1])
                    shomate['E'] = float(sline[sline.index('E')+1])
                    shomate['F'] = float(sline[sline.index('F')+1])
                    shomate['G'] = float(sline[sline.index('G')+1])
                    shomate['H'] = float(sline[sline.index('H')+1])
                    shomate['dH'] = float(sline[sline.index('dH')+1])
                    shomate['dS'] = float(sline[sline.index('dS')+1])
                    thermo['data'] = shomate
#                elif thermotype == 'standardstate':
                elif thermotype == 'NASApoly':
                    nasapoly = {}
                    sline = sl[2].split()
                    nasapoly['a1'] = float(sline[sline.index('a1')+1])
                    nasapoly['a2'] = float(sline[sline.index('a2')+1])
                    nasapoly['a3'] = float(sline[sline.index('a3')+1])
                    nasapoly['a4'] = float(sline[sline.index('a4')+1])
                    nasapoly['a5'] = float(sline[sline.index('a5')+1])
                    nasapoly['a6'] = float(sline[sline.index('a6')+1])
                    nasapoly['a7'] = float(sline[sline.index('a7')+1])
                    thermo['data'] = nasapoly
                SpeciesList[idx].thermo = thermo
                       
            if flag_bindinfo == 1:
            # extract species thermodynamic data
                if line.split() == []:
                    continue
                sl = line.split(s_seperate)
                # name of species // type of data // thermodynamic data
                spe = sl[0].replace(' ', '')
                idx = MF_MKM.get_index_species(spe, SpeciesList)
                bindsite = sl[1].replace(' ', '')
                bindenergy = float(sl[2])
                bindinfo = {}
                bindinfo['BindSite'] = bindsite
                bindinfo['BE'] = bindenergy
                SpeciesList[idx].BindInfo = bindinfo
               
            if flag_reactionthermo == 1:
            # extract reaction thermodynamic data
                if line.split() == []:
                    continue
                sl = line.split(s_seperate)
                rxn = sl[0].replace(' ', '')
                idx = MF_MKM.get_index_reactions(rxn, ReactionList)
                thermotype = sl[1].replace(' ', '')
                kinline = sl[2].split()
                if thermotype == 'Arr':
                    param = dict()
                    param['A'] = float(kinline[kinline.index('A')+1])
                    param['Ea'] = float(kinline[kinline.index('Ea')+1])
                    param['n'] = float(kinline[kinline.index('n')+1])
                elif thermotype == 'Shomate':
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
                kinetic = {}
                kinetic['type'] = thermotype
                kinetic['data'] = param
                ReactionList[idx].kinetic = kinetic
                               
            if flag_condition == 1:
            # extract reaction condition 
                if line.split() == []:
                    continue
                sl = line.split()
                if sl[0] == 'Num_condition':
                    num_condition = int(sl[1])
                else:
                    if flag_inlet == 1:
                        # name of condition // Temperature // TotalPressure  //  Total Flowrate (s-1) // Total active site (mol) // Simulation time (s) // Name_of_species  partial_pressure
                        condition = MF_MKM.ReactionCondition()
                        sl = line.split(s_seperate)
                        condi = sl[0].replace(' ', '')
                        Tem = float(sl[1])
                        T_pressure = float(sl[2])
                        T_flowrate = float(sl[3])
                        T_active = float(sl[4])
                        T_simulationtime = float(sl[5])
                        PartialP_dict = {}
                        for spe in SpeciesList:
                            if spe.phase == 'gaseous':
                                PartialP_dict[str(spe)] = 0
                        for ss in sl[6:]:
                            spe = ss.split()[0]
                            press = float(ss.split()[1])
                            PartialP_dict[spe] = press
                        condition.name = condi
                        condition.Temperature = Tem
                        condition.TotalPressure = T_pressure
                        condition.TotalFlow = T_flowrate
                        condition.TotalSite = T_active
                        condition.PartialPressure = PartialP_dict
                        condition.SimulationTime = T_simulationtime
                        
                        ConditionList.append(condition)
                    if flag_outlet == 1:
                    # name of condition // Name _of_species  rate (s-1)
                        sl = line.split(s_seperate)
                        condi = sl[0].replace(' ', '')
                        TOF_dict = {}
                        for ss in sl[1:]:
                            spe = ss.split()[0]
                            rate = float(ss.split()[1])
                            TOF_dict[spe] = rate
                        idx = MF_MKM.get_index_condition(condi, ConditionList)
                        ConditionList[idx].TurnOverFrequency = TOF_dict
                        
            if flag_parameter == 1:
                if flag_speciesparameter == 1:
                # extract parameter of the 
                    if line.split() == []:
                        continue
                    sl = line.split(s_seperate)
                    name = sl[0].replace(' ', '')
                    lb = float(sl[1])
                    ub = float(sl[2])
                    bound = [lb, ub]
                    param = {}
                    param['name'] = name
                    param['bound'] = bound
                    param['index'] = MF_MKM.get_index_species(name, SpeciesList)
                    SpeciesParameter.append(param)
                    
                if flag_reactionparameter == 1:
                # extract parameter of the 
                    if line.split() == []:
                        continue
                    sl = line.split(s_seperate)
                    name = sl[0].replace(' ', '')
                    lb = float(sl[1])
                    ub = float(sl[2])
                    bound = [lb, ub]
                    param = {}
                    param['name'] = name
                    param['bound'] = bound
                    param['index'] = MF_MKM.get_index_reactions(name, ReactionList)
                    ReactionParameter.append(param)                             
                
            
            if line.startswith(s_start):
                sl = line.split()
                section = sl[1]
                # change section flag
                if section == 'activesite':
                    flag_activesite = 1
                elif section == 'species':
                    flag_species = 1
                elif section == 'reaction':
                    flag_reaction = 1
                elif section == 'reactor':
                    flag_reactor = 1
                elif section == 'condition':
                    flag_condition = 1
                elif section == 'vibration':
                    flag_vibration = 1
                elif section == 'reaction_DFT':
                    flag_rxnDFT = 1
                elif section == 'reaction_thermo':
                    flag_reactionthermo = 1
                elif section == 'species_thermo':
                    flag_speciesthermo = 1
                elif section == 'inlet':
                    flag_inlet = 1
                elif section == 'outlet':
                    flag_outlet = 1
                elif section == 'bindinfo':
                    flag_bindinfo = 1
                elif section == 'parameter':
                    flag_parameter = 1
                elif section == 'species_parameter':
                    flag_speciesparameter = 1
                elif section == 'reaction_parameter':
                    flag_reactionparameter = 1
                else:
                    print 'Line %d: Cannot Find Section!!' %count

        
    fp.close()
    

    return ActiveSite, SpeciesList, ReactionList, ConditionList, SpeciesParameter, ReactionParameter
            
    


def is_number(s):
#    print(s)
    try:
        int(s)
        return True
    except ValueError:
#        raise
        return False
    
if __name__ == '__main__':
    a, b, c = read_file('INPUTWGS')
    for jj in b:
        print jj, jj
    for jj in c:
        print jj, jj.prefactor
                           


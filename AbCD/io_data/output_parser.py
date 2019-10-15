from AbCD import ActiveSite, GasSpecies, SurfaceSpecies, Reaction
from AbCD import CSTRCondition, VacTPDcondition, BATCHcondition
from AbCD.utils import get_index_site, get_index_species, get_index_reaction
import numpy as np

class Out_data(object):

    @staticmethod
    def print_network(network, Tem=298.15, prints=True):
        out = ''
        
        specieslist = network.specieslist
        
        # Species Basic Info
        out += 'Species Basic ThermoInfo @ %.2f K \n' % Tem
        out += '{0:^6s} {1:^10s} {2:^10s} {3:^6s} {4:^10s} {5:^12s} {6:^12s} \n'\
                .format('#', 'species', 'phase', 'mass', 'denticity', 'H(kJ/mol)', 'S(J/mol-K)')
        out += ('=='*40 + '\n')
        for i, spe in enumerate(specieslist):
            H = spe.Enthalpy(Tem)
            S = spe.Entropy(Tem)
            if spe.phase == 'gaseous':
                den = '--'
            else:
                den = str(spe.denticity)
            out += '{0:^6d} {1:^10s} {2:^10s} {3:^6.1f} {4:^10s} {5:^12.2f} {6:^12.2f} \n'\
                    .format(i, str(spe), spe.phase, spe.mass, den, H, S)

        out += '\n\n'

        # Reaction Basic Info
        out += 'Reaction Basic ThermoInfo @ %.2f K \n' % Tem
        out += '{0:^6s} {1:^6s} {2:^30s} {3:^16s} {4:^12s} {5:^6s}\n'\
                .format('#', 'name', 'reaction', 'pref (s-1)', 'Ea (kJ/mol)', 'n')
        out += ('=='*40 + '\n')
        reactionlist = network.reactionlist
        for j, rxn in enumerate(reactionlist):
            Arr = rxn.Arrhenius(Tem)
            Ea = Arr['Ea']
            n = Arr['n']
            A = Arr['A']
            out += '{0:^6d} {1:^6s} {2:^30s} {3:^16.3e} {4:^12.2f} {5:^6.2f}\n'.format(j, rxn.name, str(rxn), A, Ea, n)
        if prints:
            print(out)
        return out

    @staticmethod
    def print_optimization(network, dE_opt, dE_prime=None, Tem=298.15):
        dEa_index = network.dEa_index
        dBE_index = network.dBE_index
        Ea = dE_opt[:network._NEa]
        BE = dE_opt[network._NEa:]
#        print(Ea, BE)
#        print(dEa_index, dBE_index)
        if dE_prime is not None:
            Ea_prime = dE_prime[:network._NEa]
            BE_prime = dE_prime[network._NEa:]
        specieslist = network.specieslist
        
        out = ''
        # Species Basic Info
        out += 'Species ThermoInfo @ %.2f K (OPTIMIZED)\n' % Tem
        out += '{0:^6s} {1:^12s} {2:^12s} {3:^12s} {4:^12s} \n'\
                .format('#', 'species', 'H(kJ/mol)', 'Opt(kJ/mol)', 'Mean(kJ/mol)')
        out += ('=='*20 + '\n')
        for i, spe in enumerate(specieslist):
            H = spe.Enthalpy(Tem)
            if i not in dBE_index:
                opt = '--'
                mean = '--'
                spe.dE = 0
            else:
                opt = '%.2f' %BE[dBE_index.index(i)]
                spe.dE = BE[dBE_index.index(i)]
                if dE_prime is not None:
                    mean = '%.2f' %BE_prime[dBE_index.index(i)]
                else:
                    mean = '--'
            out += '{0:^6d} {1:^12s} {2:^12.1f} {3:^12s} {4:^12s} \n'\
                    .format(i, str(spe), H, opt, mean)

        out += '\n\n'
        # Reaction Basic Kinetic Information
        out += 'Reaction Basic Kinetic @ %.2f K (OPTIMIZED)\n' % Tem
        out += '{0:^6s} {1:^6s} {2:^30s} {3:^12s} {4:^12s} {5:^12s}\n'\
                .format('#', 'name', 'reaction', 'Ea (kJ/mol)', 'Opt(kJ/mol)',  'Mean(kJ/mol)')
        out += ('=='*40 + '\n')
        reactionlist = network.reactionlist
        for j, rxn in enumerate(reactionlist):
            Arr = rxn.Arrhenius(Tem)
            Ea_r = Arr['Ea']
            if j not in dEa_index:
                opt = '--'
                mean = '--'
            else:
                opt = '%.2f' %Ea[dEa_index.index(j)]
                if dE_prime is not None:
                    mean = '%.2f' %Ea_prime[dEa_index.index(j)]
                else:
                    mean = '--'
            out += '{0:^6d} {1:^6s} {2:^30s} {3:^12.2f} {4:^12s} {5:^12s}\n'.format(j, rxn.name, str(rxn), Ea_r, opt, mean)

            
            
        out += '\n\n'
        # Reaction Thermodynamic Information
        out += 'Reaction Basic ThermoInfo @ %.2f K (OPTIMIZED)\n' % Tem
        out += '{0:^6s} {1:^6s} {2:^30s} {3:^12s} {4:^12s} {5:^12s}\n'\
                .format('#', 'name', 'reaction', 'H0 (kJ/mol)', 'dH(kJ/mol)', 'dS(J/mol)')
        out += ('=='*40 + '\n')
        reactionlist = network.reactionlist
        for j, rxn in enumerate(reactionlist):
            
            dH = rxn.dH_Enthalpy(Tem)
            dS = rxn.dS_Entropy(Tem)
            dG0 = dH - Tem * dS/1000
            ddH = rxn.deltaH()
            
            dG_ = dG0 + ddH
            out += '{0:^6d} {1:^6s} {2:^30s} {3:^12.2f} {4:^12.2f} {5:^12.2f}\n'.format(j, rxn.name, str(rxn), dH, ddH, dS)

        print(out)

    @staticmethod
    def check_reaction(network, reaction, dE_opt=[], Tem=298.15):
        dEa_index = network.dEa_index
        dBE_index = network.dBE_index
        Ea = dE_opt[:network._NEa]
        BE = dE_opt[network._NEa:]
        specieslist = network.specieslist
        
        out = ''
        for i, spe in enumerate(specieslist):
            H = spe.Enthalpy(Tem)
            if i not in dBE_index:
                spe.dE = 0
            else:
                spe.dE = BE[dBE_index.index(i)]
        for reactant in reaction.reactant:
            react = reactant[0]
            stoi = reactant[1]
            out += '{0:^12s} {1:^6.1f} {2:^12.2f} {3:^8.2f} {4:^6.2f}\n'.format(react, stoi, react.Enthalpy(Tem), react.Entropy(Tem), react.dE)
        for reactant in reaction.product:
            react = reactant[0]
            stoi = reactant[1]
            out += '{0:^12s} {1:^6.1f} {2:^12.2f} {3:^8.2f} {4:^6.2f}\n'.format(react, stoi, react.Enthalpy(Tem), react.Entropy(Tem), react.dE)
        
        dH0 = reaction.dH_Enthalpy(Tem)
        dS0 = reaction.dS_Entropy(Tem)
        
        dG0 = dH0 - Tem * dS0/1000
        ddH = reaction.deltaH()
        
        dHf = dH0 + ddH
        dGf = dHf - Tem * dS0/1000
        import math
        K0 = math.exp( -dG0 * 1000/ (Tem * 8.314))
        Kf = math.exp( -dGf * 1000/ (Tem * 8.314))
        
        out += 'Delta H = %.2f kJ/mol\n' %(ddH)
        out += '{0:^6s} {1:^30s} {2:^12.2f} {3:^12.2f} {4:^12.2f} {5:^12.4e}\n'.\
                format(reaction.name, str(reaction), dH0, dS0, dG0, K0)
        out += '{0:^6s} {1:^30s} {2:^12.2f} {3:^12.2f} {4:^12.2f} {5:^12.4e}\n'.\
                format('--', '--', dHf, dS0, dGf, Kf)
        
        print(out)

    @staticmethod
    def toVizPES(network, dE, barListIn=[], connListIn=[], 
                 mainBarIndex=[], sideBarIndex=[],
                 mainConnIndex=[], sideConnIndex=[],
                 neglect_h='False',
                 Tem=298.15, e_type='H', shift=True,
                 pesname='pesFile.pes'):
        '''
        Write network information to CataViz for Potential Energy Surface 
        '''

        assign_correction(network, dE)
        

        # CREATE BARLIST
        barList = []
        E0 = 0
        for i, indices_dict in enumerate(barListIn):
            E = 0
            name = []
            
            for idx in indices_dict.keys():
                stoi = indices_dict[idx]
                spe = network.specieslist[idx]
                E += stoi * spe.Enthalpy(Tem, corr=True)
                if not (neglect_h and str(spe) == 'H*'):
                    str_stoi = str(stoi) if stoi != 1 else ''
                    name.append(str_stoi + spe.unicode_repr())
            if i == 0 and shift:
                E0 = E
            E = E - E0
            name = '\n+'.join(name)
            barList.append({'E': E, 'name':name})
        
        connectorList = []
        for i, indices in enumerate(connListIn):
            name = []
            
            # Needed Setting
            idx = indices['index']
            left = indices['left']
            right = indices['right']
            dirc = indices.get('dirc', 1)
            # Default Setting
            shape = indices.get('shape', 'quad')
            
            rxn = network.reactionlist[idx]
            E = barList[left]['E'] + rxn.deltaEa(Tem, dirc)
            
            left = barList[left]['name']
            right = barList[right]['name']
            
            name = rxn.name + '=' + left + '_' + 'right'
            connectorList.append({'E': E, 'name':name, 'left': left, 'right': right, 'shape': shape})

        mainBarName = [barList[i]['name'] for i in mainBarIndex]
        sideBarName = [barList[i]['name'] for i in sideBarIndex]            
        mainConnName = [connectorList[i]['name'] for i in mainConnIndex]
        sideConnName = [connectorList[i]['name'] for i in sideConnIndex]
        
        pes = {}
        pes["barList"] = barList
        pes["connectList"] = connectorList
        pes["mainBar"] = mainBarName
        pes["mainConnect"] = mainConnName
        pes["sideConnect"] = sideConnName
        
        # TO Json file
        import json
        with open(pesname, 'w') as fp:
            json.dump(pes, fp, indent=4)

    @staticmethod
    def toVizPES_percentile(network, dEs, percentile, barListIn=[], connListIn=[], 
                 mainBarIndex=[], sideBarIndex=[],
                 mainConnIndex=[], sideConnIndex=[],
                 setting={}, N=1000,
                 Tem=298.15, e_type='H', shift=True,
                 pesname='pesFile.pes'):
        '''
        Write network information to CataViz for Potential Energy Surface 
        '''
        assert percentile >= 0 and percentile <= 100
        ndata, _ = dEs.shape
        dn = int(ndata / N)
        bar_ensemble = []; connect_ensemble = []
        for dE in dEs[::dn, :]:
            assign_correction(network, dE)

            # CREATE BARLIST
            barE, barList = [], []
            E0 = 0
            for i, indices_dict in enumerate(barListIn):
                E = 0
                for idx in indices_dict.keys():
                    stoi = indices_dict[idx]
                    spe = network.specieslist[idx]
                    E += stoi * spe.Enthalpy(Tem, corr=True)
                if i == 0 and shift:
                    E0 = E
                E = E - E0
                barE.append(E)
                barList.append({'E': E, 'name': ""})
            bar_ensemble.append(barE)
            
            connectE = []
            for i, indices in enumerate(connListIn):
                # Needed Setting
                idx = indices['index']
                left = indices['left']
                right = indices['right']
                dirc = indices.get('dirc', 1)
                # Default Setting
                shape = indices.get('shape', 'quad')
                
                rxn = network.reactionlist[idx]
                E = barList[left]['E'] + rxn.deltaEa(Tem, dirc)
                connectE.append(E)
            connect_ensemble.append(connectE)

        # get the percentile energy
        bar_ensemble = np.array(bar_ensemble)
        connect_ensemble = np.array(connect_ensemble)
        bar_percent = np.percentile(bar_ensemble, percentile, axis=0)
        connect_percent = np.percentile(connect_ensemble, percentile, axis=0)
        # build output

        
        # CREATE BARLIST
        barList = []
        E0 = 0
        for i, indices_dict in enumerate(barListIn):
            E = 0
            name = []
            for idx in indices_dict.keys():
                stoi = indices_dict[idx]
                spe = network.specieslist[idx]
                str_stoi = str(stoi) if stoi != 1 else ''
                name.append(str_stoi + spe.name)
            E = bar_percent[i]
            name = '+\n'.join(name)
            barList.append({'E': E, 'name':name})
        
        connectorList = []
        for i, indices in enumerate(connListIn):
            name = []
            
            # Needed Setting
            idx = indices['index']
            left = indices['left']
            right = indices['right']
            dirc = indices.get('dirc', 1)
            # Default Setting
            shape = indices.get('shape', 'quad')
            rxn = network.reactionlist[idx]
            E = connect_percent[i]
            left = barList[left]['name']
            right = barList[right]['name']
            
            name = rxn.name + '=' + left + '_' + 'right'
            connectorList.append({'E': E, 'name':name, 'left': left, 'right': right, 'shape': shape})

        mainBarName = [barList[i]['name'] for i in mainBarIndex]
        sideBarName = [barList[i]['name'] for i in sideBarIndex]            
        mainConnName = [connectorList[i]['name'] for i in mainConnIndex]
        sideConnName = [connectorList[i]['name'] for i in sideConnIndex]
        
        pes = {}
        pes["barList"] = barList
        pes["connectList"] = connectorList
        pes["mainBar"] = mainBarName
        pes["mainConnect"] = mainConnName
        pes["sideConnect"] = sideConnName
        
        # TO Json file
        import json
        with open(pesname, 'w') as fp:
            json.dump(pes, fp, indent=4)


def assign_correction(network, dE):
    dEa_index = network.dEa_index
    dBE_index = network.dBE_index
    Ea = dE[:network._NEa]
    BE = dE[network._NEa:]
    # Assign correction on species
    for i, spe in enumerate(network.specieslist):
        if i not in dBE_index:
            spe.dE = 0
        else:
            spe.dE = BE[dBE_index.index(i)]

    # Assign correction on reactions
    for j, rxn in enumerate(network.reactionlist):
        if j not in dEa_index:
            rxn.dE = 0
        else:
            rxn.dE = Ea[dEa_index.index(j)]
from AbCD import ActiveSite, GasSpecies, SurfaceSpecies, Reaction
from AbCD import CSTRCondition, VacTPDcondition, BATCHcondition
from AbCD.utils import get_index_site, get_index_species, get_index_reaction


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
                    .format(i+1, str(spe), spe.phase, spe.mass, den, H, S)

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
            out += '{0:^6d} {1:^6s} {2:^30s} {3:^16.3e} {4:^12.2f} {5:^6.2f}\n'.format(j+1, rxn.name, str(rxn), A, Ea, n)
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
                    .format(i+1, str(spe), H, opt, mean)

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
            out += '{0:^6d} {1:^6s} {2:^30s} {3:^12.2f} {4:^12s} {5:^12s}\n'.format(j+1, rxn.name, str(rxn), Ea_r, opt, mean)

            
            
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
            out += '{0:^6d} {1:^6s} {2:^30s} {3:^12.2f} {4:^12.2f} {5:^12.2f}\n'.format(j+1, rxn.name, str(rxn), dH, ddH, dS)

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
    def toCataViz(network, dE, main_bars_index=[], main_connectors_index=[],
                  side_bars_index=[], setting={}):
        '''
        Write network information to CataViz for Potential Energy Surface 
        '''
        pass
        
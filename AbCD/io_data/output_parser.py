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
            else:
                opt = '%.2f' %BE[dBE_index.index(i)]
                if dE_prime is not None:
                    mean = '%.2f' %BE_prime[dBE_index.index(i)]
                else:
                    mean = '--'
            out += '{0:^6d} {1:^12s} {2:^12.1f} {3:^12s} {4:^12s} \n'\
                    .format(i+1, str(spe), H, opt, mean)

        out += '\n\n'
        # Reaction Basic Info
        out += 'Reaction Basic ThermoInfo @ %.2f K (OPTIMIZED)\n' % Tem
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

        print(out)
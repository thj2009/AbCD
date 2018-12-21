# IO data for cstr model

def cstr_output(cstr):
    temp = cstr._condition.Temperature

    out = ''
    # Reaction Arrhenius Constant
    out += 'r = A*T^n * exp(-Ea/RT)  @ %.2f K' %(temp) + '\n'
    out += '{0:25s} {1:^12s} {2:^12s} {3:>12s} {4:>12s}'\
            .format('Reaction', 'Prefactor', 'n', 'oldEa', 'newEa') + '\n'
    out += '-'*80 + '\n'
    for idx, rxn in enumerate(cstr.reaction_net.reactionlist):
        out += '{0:25s} {1:^-12,.2e} {2:^-12,.2f} {3:>12,.2f} {4:>12,.2f}'\
                .format(rxn, rxn.Arrhenius(temp)['A'], rxn.Arrhenius(temp)['n'],
                        rxn.Arrhenius(temp)['Ea'], cstr.Energy['activation'][idx]) + '\n'
    out += '\n'

    # Reaction Enthalpy
    out += 'Reaction Enthalpy and Entropy @ %.2f K' %temp + '\n'
    out += '{0:25s} {1:^12s} {2:^12s} {3:^12s}'\
            .format('Reaction', 'dS(J/mol/K)', 'old_dH', 'new_dH') + '\n'
    out += '-'*80 + '\n'
    for idx, rxn in enumerate(cstr.reaction_net.reactionlist):
        out += '{0:25s} {1:^-12,.2f} {2:^-12,.2f} {3:^-12,.2f}'\
                .format(rxn, rxn.dS_Entropy(temp), rxn.dH_Enthalpy(temp),
                        cstr.Energy['enthalpy'][idx]) + '\n'
    out += '\n'

    # Species Enthalpy and Entropy
    out += 'Species Enthalpy and Entropy' + '\n'
    out += '-'*80 + '\n'
    out += '{0:25s} {1:^12s} {2:^12s}'\
            .format('Species', 'S(J/mol/K)', 'H(kJ/mol)') + '\n'
    out += '-'*80 + '\n'
    for idx, spe in enumerate(cstr.reaction_net.specieslist):
        out += '{0:25s} {1:^-12,.2f} {2:^-12,.2f}'\
                .format(spe, spe.Entropy(temp), spe.Enthalpy(temp)) + '\n'
    out += '\n'

    # Reaction Rate Partition
    out += 'Reaction Rate (s-1)' + '\n'
    out += '-' * 80 + '\n'
    out += '{0:25s} {1:^12s} {2:^12s} {3:^12s}'\
            .format('Reaction', 'Net', 'Rfor', 'Rrev') + '\n'
    out += '-' * 80 + '\n'
    for idx, rxn in enumerate(cstr.reaction_net.reactionlist):
        out += '{0:25s} {1:^-12,.2e} {2:^-12,.2e} {3:^-12,.2e}'\
                .format(rxn, cstr.Rate['rnet'][idx], cstr.Rate['rfor'][idx],
                        cstr.Rate['rrev'][idx]) + '\n'
    out += '\n'

    # Surface Coverage Parittion
    out += 'Surface Coverage (fraction)' + '\n'
    out += '-' * 80 + '\n'
    for idx, spe in enumerate(cstr.reaction_net.specieslist):
        if idx >= cstr.reaction_net.Ngas:
            out += '{0:20s} {1:^-8,.2e} '\
                    .format(spe, cstr.coverage[idx - cstr.reaction_net.Ngas]) + '\n'
    out += '{0:20s} {1:^-8,.2e}' .format('Total', sum(cstr.coverage)) + '\n'
    out += '\n'

    # Reaction rate Constant and Equilibrium Constant
    out += 'Kinetic and Equilibrium Constant' + '\n'
    out += '{0:25s} {1:^12s} {2:^12s} {3:^12s} {4:^12s}'\
            .format('Reaction', 'Keq', 'Qeq', 'kf', 'kr') + '\n'
    out += '-' * 80 + '\n'
    for idx, rxn in enumerate(cstr.reaction_net.reactionlist):
        out += '{0:25s} {1:^-12,.2e} {2:^-12,.2e} {3:>-12,.2e} {4:>-12,.2e}'\
                .format(rxn, cstr.EqCons_RCons['Keq'][idx], cstr.EqCons_RCons['Qeq'][idx],
                        cstr.EqCons_RCons['kf'][idx], cstr.EqCons_RCons['kr'][idx]) + '\n'
    out += '\n'

    # Species Enthalpy and Entropy
    out += 'Reaction Energy, Rate Constant, Equilibrium' + '\n'
    out += '-' * 120 + '\n'
    out += '{0:25s} {1:^12s} {2:^12s} {3:^12s} {4:^12s} {5:^12s} {6:^12s}'\
            .format('Reaction', 'dH(kJ/mol)', 'dS(J/mol/K)', 'dG(kJ/mol)',
                    'K', 'prefactor', 'Ea(kJ/mol)') + '\n'
    out += '-' * 120 + '\n'
    for idx, rxn in enumerate(cstr.reaction_net.reactionlist):
        dH = cstr.Energy['enthalpy'][idx]
        dS = rxn.dS_Entropy(temp)
        dG = rxn.dG_GibbsFreeEnergy(temp)
        K = cstr.EqCons_RCons['Keq'][idx]
        A = rxn.Arrhenius(temp)['A']
        Ea = cstr.Energy['activation'][idx]
        out += '{0:25s} {1:^-12,.2f} {2:^-12,.2f} {3:^-12,.2f} {4:^-12,.2e} {5:^-12,.2e} {6:^-12,.2f}'\
                .format(rxn, dH, dS, dG, K, A, Ea) + '\n'
    out += '\n'
    return out
import os
import json
from AbCD import ActiveSite, GasSpecies, SurfaceSpecies, Reaction
from AbCD import CSTRCondition, VacTPDcondition
from AbCD.utils import get_index_site, get_index_species, get_index_reaction

class In_data(object):
    '''
    Read Input data
    '''
    @staticmethod
    def read_site(site_file, _dir_):
        '''
        Read the active site information
        :param site_string: json readable string
        '''
        with open(site_file, 'r') as fp:
            datalist = json.load(fp)
        for data in datalist:
            site = ActiveSite()
            site.name = data['name']
            site.metal = data['metal']
            site.struct = data['struct']
            site.facet = data['facet']
            site.site = data['site']
            site.lattice_constant = data['lattice_constant']
            site.mass = float(data['mass'])
            _dir_['sitelist'].append(site)

    @staticmethod
    def read_species(species_file, _dir_):
        with open(species_file, 'r') as f:
            datalist = json.load(f)
        for data in datalist:
            phase = data['phase']
            if phase in ['gaseous', 'gas', 'g']:
                spe = GasSpecies('')
                spe.name = data['name']
                spe.formula = data['formula']
                spe.mass = float(data['mass'])
            elif phase in ['surface', 'surf', 's']:
                spe = SurfaceSpecies()
                spe.name = data['name']
                spe.formula = data['formula']
                spe.mass = float(data['mass'])
                idx = get_index_site(data['site_name'], _dir_['sitelist'])
                spe.site = _dir_['sitelist'][idx]
                spe.denticity = int(data['denticity'])
            _dir_['specieslist'].append(spe)

    @staticmethod
    def read_reaction(reaction_file, _dir_):
        '''
        Read the reaction mechanism given list of reactions
        :param reaction_string: json readable string
        '''
        with open(reaction_file, 'r') as f:
            datalist = json.load(f)
        for data in datalist:
            rxn = Reaction()
            rxn.name = data['name']
            reactant = data['reactant']
            product = data['product']
            cons_react = []
            for react, stoi in reactant:
                idx = get_index_species(react, _dir_['specieslist'])
                cons_react.append((_dir_['specieslist'][idx], stoi))
            cons_prod = []
            for react, stoi in product:
                idx = get_index_species(react, _dir_['specieslist'])
                cons_prod.append((_dir_['specieslist'][idx], stoi))
            rxn.reactant = cons_react
            rxn.product = cons_prod
            _dir_['reactionlist'].append(rxn)

    @staticmethod
    def read_species_thermo(species_thermo_file, _dir_):
        '''
        Assign the species thermo calculation result to species list,
        :param species_thermo_string: json readable string for list of dictionary
        :param specieslist: list stored all species information
        '''
        with open(species_thermo_file, 'r') as f:
            datalist = json.load(f)
        for data in datalist:
            spe = data['name']
            idx = get_index_species(spe, _dir_['specieslist'])
            _dir_['specieslist'][idx].thermo = data['thermo']

    @staticmethod
    def read_reaction_kinetic(reaction_kinetic_file, _dir_):
        '''
        Assign the reaction kinetic calculation result to reaction list,
        :param reaction_kinetic_string: json readable string for list of dictionary
        :param reactionlist: list stored all reaction information
        '''
        with open(reaction_kinetic_file, 'r') as f:
            datalist = json.load(f)
        for data in datalist:
            rxn = data['name']
            idx = get_index_reaction(rxn, _dir_['reactionlist'])
            _dir_['reactionlist'][idx].kinetic = data['kinetic']

    @staticmethod
    def read_reaction_dft(reaction_dft_file, _dir_):
        '''
        Assign the reaction dft calculation result to reaction list,
        :param reaction_dft_string: json readable string for list of dictionary
        :param reactionlist: list stored all reaction information
        '''
        with open(reaction_dft_file, 'r') as f:
            datalist = json.load(f)
        for data in datalist:
            rxn = data['name']
            idx = get_index_reaction(rxn, _dir_['reactionlist'])
            _dir_['reactionlist'][idx].dft_data = data['dft_data']

    @staticmethod
    def read_species_vib(species_vib_file, _dir_):
        '''
        Assign the vibrational frequency to species list,
        :param vib_string: json readable string for list of dictionary
        :param specieslist: specieslist stored all species information
        '''
        with open(species_vib_file, 'r') as f:
            datalist = json.load(f)
        for data in datalist:
            spe = data['name']
            idx = get_index_species(spe, _dir_['specieslist'])
            _dir_['specieslist'][idx].vibfreq = data['vib']
    
    @staticmethod
    def read_reaction_condition(reaction_condition_file, conditionlist):
        with open(reaction_condition_file, 'r') as f:
            datalist = json.load(f)
        for data in datalist:
            if data['type'] == 'CSTR':
                condi = CSTRCondition()
                condi.name = data['name']
                condi.Temperature = data['Temperature']
                condi.TotalPressure = data['TotalPressure']
                condi.TotalFlow = data['TotalFlow']
                condi.TotalSite = data['TotalSite']
                condi.SimulationTime = data['SimulationTime']
                condi.PartialPressure = data['PartialPressure']
                condi.TurnOverFrequency = data['TurnOverFrequency']
            elif data['type'] == 'Batch':
                pass
            elif data['type'] == 'VacTPD':
                condi = VacTPDcondition()
                condi.name = data['name']
                condi.T0 = data['T0']
                condi.Tf = data['Tf']
                condi.Beta = data['Beta']
                condi.Ntime = data['Ntime']
                condi.InitCoverage = data['InitCoverage']
                if 'SimulationTime' not in data.keys():
                    condi._calSim()
                else:
                    condi.SimulationTime = data['SimulationTime']
                if 'TemGrid' not in data.keys():
                    condi._calGrid()
                else:
                    condi.TimeGrid = data['TimeGrid']
                    condi.TemGrid = data['TemGrid']
                if 'TemProfile' in data.keys():
                    condi.TemProfile = data['TemProfile']
                    condi.RateProfile = data['RateProfile']
            conditionlist.append(condi)
        return conditionlist

    @staticmethod
    def read_parameter():
        # TODO
        pass

    @staticmethod
    def load_condition(data_fdr):
        conditionlist = []
        method_to_call = getattr(In_data, 'read_reaction_condition')
        method_to_call(os.path.join(data_fdr, 'reaction_condition.mkm'),
                       conditionlist)
        return conditionlist

    @staticmethod
    def load_mkm(data_fdr):
        '''
        Read composite mkmfile
        return: active site information, species info, reaction info, and full reaction
        network
        '''
        default = [
            'activesite.mkm',
            'species.mkm',
            'reaction.mkm',
            'species_thermo.mkm',
            'reaction_kinetic.mkm',
            'species_vibration.mkm',
            'reaction_dft.mkm'
            ]
        functions = ['read_site',
                     'read_species',
                     'read_reaction',
                     'read_species_thermo',
                     'read_reaction_kinetic',
                     'read_species_vib',
                     'read_reaction_dft'
                     ]
        _dir_ = {'sitelist': [], 'specieslist': [], 'reactionlist': []}
        for f, attr in zip(default, functions):
            method_to_call = getattr(In_data, attr)
            method_to_call(os.path.join(data_fdr, f), _dir_)
        return _dir_['sitelist'], _dir_['specieslist'], _dir_['reactionlist']


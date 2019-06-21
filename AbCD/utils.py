# Utility function

class Constant:

    c = 299792458               # Speed of light m s-1
    kb = 1.38064852e-23         # Boltzmann Constant J K-1
    h = 6.62607004e-34          # Blankc Constant   J s
    NA = 6.02214086e23          # Avogadro Constant mol-1
    Rg = 8.3144598              # Ideal Gas Constant

    # Converter
    hatree2kJmol = 2625.5       # Hatree to KJ.mol
    hatree2J = 4.35974e-18      # Hatree to J
    half_hc = 1/2. * h * c * 100 * NA/1000. # kJ cm /mol
    u_2_kg = 1.66054e-27

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def get_index_species(spe, spelist, neglect=False):
    for i in range(len(spelist)):
            if spe == str(spelist[i]):
                return i
    if not neglect:
        raise NameError(spe + ' is not in the species list')


def get_index_reaction(rxnname, rxnlist):
    for i in range(len(rxnlist)):
        if rxnname == rxnlist[i].name:
            return i
    raise NameError(rxnname+ ' is not in the reaction list')

def get_index_site(sitename, sitelist):
    for i in range(len(sitelist)):
        if sitename == sitelist[i].name:
            return i
    raise NameError(sitename+ ' is not in the site list')

def get_index_condition(condiname, condilist):
    for i in range(len(condilist)):
        if condiname == condilist[i].name:
            return i
    raise NameError(condiname+ ' is not in the condition list')

def check_index(index, checklist):
    for i in checklist:
        if i == index:
            return True
    return False

def find_index(index, checklist):
    for i in range(len(checklist)):
        if index == checklist[i]:
            return i

def other_reaction(specieslist, stoi, Tem=298.75):
    H = 0
    S = 0
    for idx, stoi_ in stoi:
        H += reactionlist[idx].Enthalpy(Tem) * stoi_
        S += reactionlist[idx].Entropy(Tem) * stoi_
    return H, S
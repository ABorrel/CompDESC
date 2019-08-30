from rdkit.Chem import Descriptors

def getFpDensityMorgan1(mol):
    return Descriptors.FpDensityMorgan1(mol)

def getFpDensityMorgan2(mol):
    return Descriptors.FpDensityMorgan2(mol)

def getFpDensityMorgan3(mol):
    return Descriptors.FpDensityMorgan3(mol)


_morgan = {"FpDensityMorgan1": getFpDensityMorgan1,
           "FpDensityMorgan2": getFpDensityMorgan2,
           "FpDensityMorgan3": getFpDensityMorgan3}



def GetMorgan(mol):
    result = {}
    for DesLabel in _morgan.keys():
        result[DesLabel] = _morgan[DesLabel](mol)
    return result

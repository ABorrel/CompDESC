# Transform by Alexandre Borrel from PYDPI for python 3.6 with rdkit version 2019-3

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, Crippen

#### General functions ####
def getCountByElementNumber(mol, AtomicNumber=6):
    out = 0.0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNumber:
            out = out + 1
    return out

def getPathMolecule(mol, PathLength=2):
    return len(Chem.FindAllPathsOfLengthN(mol,PathLength,useBonds=1))

## desriptor ##
def getnH(mol):
    Hmol = Chem.AddHs(mol)
    out = getCountByElementNumber(Hmol, 1)
    return out

def getHalCount(mol):
    out=0.0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==9 or atom.GetAtomicNum()==17 or atom.GetAtomicNum()==35 or atom.GetAtomicNum()==53 or atom.GetAtomicNum()==85 or atom.GetAtomicNum()==117:
            out=out+1
    return out

def getNumHeteroatoms(mol):
    return Descriptors.NumHeteroatoms(mol)

def getHeavyAtomCount(mol):
    return mol.GetNumHeavyAtoms()

def getFCount(mol):
    return getCountByElementNumber(mol, 9)

def getClCount(mol):
    return getCountByElementNumber(mol, 17)

def getBrCount(mol):
    return getCountByElementNumber(mol, 35)

def getICount(mol):
    return getCountByElementNumber(mol, 53)

def getCCount(mol):
    return getCountByElementNumber(mol, 6)

def getPCount(mol):
    return getCountByElementNumber(mol, 15)

def getSCount(mol):
    return getCountByElementNumber(mol, 16)

def getOCount(mol):
    return getCountByElementNumber(mol, 8)

def getNCount(mol):
    return getCountByElementNumber(mol, 7)

def getRingCount(mol):
    return Descriptors.RingCount(mol)

def getNumRotatableBonds(mol):
    return Descriptors.NumRotatableBonds(mol)

def getNumHDonors(mol):
    return Descriptors.NumHDonors(mol)

def getNumHAcceptors(mol):
    return Descriptors.NumHAcceptors(mol)

def getSingleBoundCount(mol):
    out = 0.0
    for bond in mol.GetBonds():
        if bond.GetBondType().name=='SINGLE':
            out = out + 1
    return out

def getDoubleBoundCount(mol):
    out = 0.0
    for bond in mol.GetBonds():
        if bond.GetBondType().name=='DOUBLE':
            out = out + 1
    return out

def getTripleBoundCount(mol):
    out = 0.0
    for bond in mol.GetBonds():
        if bond.GetBondType().name == 'TRIPLE':
            out = out + 1
    return out

def getArBoundCount(mol):
    out = 0.0
    for bond in mol.GetBonds():
        if bond.GetBondType().name == 'AROMATIC':
            out = out + 1
    return out

def getNumAllatoms(mol):
    return Chem.AddHs(mol).GetNumAtoms()

def getPath1Count(mol):
    return getPathMolecule(mol, 1)

def getPath2Count(mol):
    return getPathMolecule(mol, 2)

def getPath3Count(mol):
    return getPathMolecule(mol, 3)

def getPath4Count(mol):
    return getPathMolecule(mol, 4)

def getPath5Count(mol):
    return getPathMolecule(mol, 5)

def getPath6Count(mol):
    return getPathMolecule(mol, 6)

def getNHOHCount(mol):
    return Descriptors.NHOHCount(mol)

def getNOCount(mol):
    return Descriptors.NOCount(mol)

def getNumAliphaticCarbocycles(mol):
    return Descriptors.NumAliphaticCarbocycles(mol)

def getNumAliphaticHeterocycles(mol):
    return Descriptors.NumAliphaticHeterocycles(mol)

def getNumAliphaticRings(mol):
    return Descriptors.NumAliphaticRings(mol)

def getNumAromaticCarbocycles(mol):
    return Descriptors.NumAromaticCarbocycles(mol)

def getNumAromaticHeterocycles(mol):
    return Descriptors.NumAromaticHeterocycles(mol)

def getNumAromaticRings(mol):
    return Descriptors.NumAromaticRings(mol)

def getNumSaturatedCarbocycles(mol):
    return Descriptors.NumSaturatedCarbocycles(mol)

def getNumSaturatedHeterocycles(mol):
    return Descriptors.NumSaturatedHeterocycles(mol)

def getNumSaturatedRings(mol):
    return Descriptors.NumSaturatedRings(mol)

def getFractionCSP3(mol):
    return Descriptors.FractionCSP3(mol)

def getLipinskiHBA(mol):
    return rdMolDescriptors.CalcNumLipinskiHBA(mol)

def getLipinskiHBD(mol):
    return rdMolDescriptors.CalcNumLipinskiHBD(mol)

def getStereoCenters(mol):
    return rdMolDescriptors.CalcNumAtomStereoCenters(mol)

def getUnspecifiedStereoCenters(mol):
    return rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)

def getAmideBonds(mol):
    return rdMolDescriptors.CalcNumAmideBonds(mol)

def getLipinskiFail(mol):
    '''
    Returns which of Lipinski's rules a molecule has failed, or an empty list
    
    Lipinski's rules are:
    Hydrogen bond donors <= 5
    Hydrogen bond acceptors <= 10
    Molecular weight < 500 daltons
    logP < 5
    '''
    passed = []
    failed = []
    
    num_hdonors = Lipinski.NumHDonors(mol)
    num_hacceptors = Lipinski.NumHAcceptors(mol)
    mol_weight = Descriptors.MolWt(mol)
    mol_logp = Crippen.MolLogP(mol)
    
    failed = []
    
    if num_hdonors > 5:
        failed.append('Over 5 H-bond donors, found %s' % num_hdonors)
    else:
        passed.append('Found %s H-bond donors' % num_hdonors)
        
    if num_hacceptors > 10:
        failed.append('Over 10 H-bond acceptors, found %s' \
        % num_hacceptors)
    else:
        passed.append('Found %s H-bond acceptors' % num_hacceptors)
        
    if mol_weight >= 500:
        failed.append('Molecular weight over 500, calculated %s'\
        % mol_weight)
        
    if mol_logp >= 5:
        failed.append('Log partition coefficient over 5, calculated %s' \
        % mol_logp)
    else:
        passed.append('Log partition coefficient: %s' % mol_logp)
    
    return len(failed)




##################
# Main function  #
##################

_constitutional={"nH": getnH,
                 "HalCount": getHalCount,
                 "NumHeteroatoms": getNumHeteroatoms,
                 "HeavyAtomCount": getHeavyAtomCount,
                 "Fcount": getFCount,
                 "ClCount": getClCount,
                 "BrCount": getBrCount,
                 "ICount": getICount,
                 "CCount": getCCount,
                 "PCount": getPCount,
                 "SCount": getSCount,
                 "OCount": getOCount,
                 "NCount": getNCount,
                 "RingCount": getRingCount,
                 "NumRotatableBonds": getNumRotatableBonds,
                 "NumHDonors": getNumHDonors,
                 "NumHAcceptors": getNumHAcceptors,
                 "SingleBoundCount": getSingleBoundCount,
                 "DoubleBoundCount": getDoubleBoundCount,
                 "ArBoundCount": getArBoundCount,
                 "TripleBoundCount": getTripleBoundCount,
                 "NumAllatoms": getNumAllatoms,
                 "Path1Count": getPath1Count,
                 "Path2Count": getPath2Count,
                 "Path3Count": getPath3Count,
                 "Path4Count": getPath4Count,
                 "Path5Count": getPath5Count,
                 "Path6Count": getPath6Count,
                 "NHOHCount": getNHOHCount,
                 "NOCount": getNOCount,
                 "NumAliphaticCarbocycles": getNumAliphaticCarbocycles,
                 "NumAliphaticHeterocycles": getNumAliphaticHeterocycles,
                 "NumAliphaticRings": getNumAliphaticRings,
                 "NumAromaticCarbocycles": getNumAromaticCarbocycles,
                 "NumAromaticHeterocycles": getNumAromaticHeterocycles,
                 "NumAromaticRings": getNumAromaticRings,
                 "NumSaturatedCarbocycles": getNumSaturatedCarbocycles,
                 "NumSaturatedHeterocycles": getNumSaturatedHeterocycles,
                 "NumSaturatedRings": getNumSaturatedRings,
                 "FractionCSP3": getFractionCSP3,
                 "NumLipinskiHBA": getLipinskiHBA,
                 "NumLipinskiHBD": getLipinskiHBD,
                 "NumStereocenters": getStereoCenters,
                 "NumUnspecifiedStereocenters": getUnspecifiedStereoCenters,
                 "NumAmideBonds":getAmideBonds,
                 "LipinskiFail":getLipinskiFail}


def GetConstitutional(mol):
    dresult={}
    for DesLabel in _constitutional.keys():
        dresult[DesLabel] = round(_constitutional[DesLabel](mol),6)
    return dresult


from rdkit import Chem
from rdkit.Chem import Descriptors

from rdkit.Chem import Crippen
from rdkit.Chem import MolSurf as MS
from rdkit.Chem import rdMolDescriptors

import math


def getMolLogP(mol):
    return Descriptors.MolLogP(mol)

def getMolLogP2(mol):
    logP = Descriptors.MolLogP(mol)
    return logP**2

def getMolMR(mol):
    return Descriptors.MolMR(mol)

def getHy(mol):
    """
    #################################################################
    Calculation of hydrophilicity factor. The hydrophilicity
    index is described in more detail on page 225 of the
    Handbook of Molecular Descriptors (Todeschini and Consonni 2000).
    In first hydrogen are added in the molecules
    => error of 0.01 compare to the original table
    #################################################################
    """
    molH = Chem.AddHs(mol)
    A = float(molH.GetNumHeavyAtoms())
    if A == 0.0:
        print("Chem Error: (molproperty l.35) No heavy atom")# see to put that in a log file
        return 0.0
    NC = 0.0
    NHy = 0.0# only H connected to N, O, S
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6:
            NC = NC + 1
        elif atom.GetAtomicNum() == 7 or atom.GetAtomicNum() == 8 or atom.GetAtomicNum() == 16:
            nHtemp = 0  # to be counted the group should have only one H
            atomn = atom.GetNeighbors()
            for n in atomn:
                if n.GetAtomicNum() == 1:
                    nHtemp = nHtemp + 1
            if nHtemp == 1:
                NHy = NHy + 1
    #print(NHy, A, NC)
    numerator = (1 + NHy) * math.log((1 + NHy), 2) + NC * ((1.0 / A) * (math.log(1.0 / A, 2))) + math.sqrt((NHy) / A**2)
    demonimator = math.log(1+A, 2)
    out = numerator / demonimator
    return out

def getBondNumber(mol,bondtype='SINGLE'):
    out = 0
    for bond in mol.GetBonds():
        if bond.GetBondType().name==bondtype:
            out = out + 1
    return out


def getUI(mol):
    nd=getBondNumber(mol,bondtype='DOUBLE')
    nt=getBondNumber(mol,bondtype='TRIPLE')
    na=getBondNumber(mol,bondtype='AROMATIC')
    res=math.log((1.0+nd+nt+na),2)
    return res

def getHeavyAtomMolWt(mol):
    return Descriptors.HeavyAtomMolWt(mol)

def getAWeightHeavyAtom(mol):
    mol = Chem.RemoveHs(mol)# remove H
    nbatom = mol.GetNumAtoms()
    molWt = Descriptors.HeavyAtomMolWt(mol)
    avgWt = molWt/nbatom
    return avgWt

def getExactMolWt(mol):
    return Descriptors.ExactMolWt(mol)

def getMolWt(mol):
    return Descriptors.MolWt(mol)

def getqed(mol):
    return Descriptors.qed(mol)


_molProperty={"MolLogP": getMolLogP,
              "MolLogP2": getMolLogP2,
              "MolMR": getMolMR,
              "Hy": getHy,
              "UI": getUI,
              "HeavyAtomMolWt": getHeavyAtomMolWt,
              "AWeightHeavyAtom": getAWeightHeavyAtom,
              "ExactMolWt": getExactMolWt,
              "MolWt": getMolWt,
              "qed": getqed}


def GetMolecularProperty(mol):
    dresult={}
    for DesLabel in _molProperty.keys():
        dresult[DesLabel] = round(_molProperty[DesLabel](mol), 6)
    return dresult

from rdkit import Chem
from rdkit.Chem import MolSurf as MOE , Descriptors
from rdkit.Chem.EState import EState_VSA as EVSA

def getLabuteASA(mol):
    return Descriptors.LabuteASA(mol)

def getTPSA(mol):
    return Descriptors.TPSA(mol)

def getSlogP_VSA1(mol):
    return Descriptors.SlogP_VSA1(mol)

def getSlogP_VSA2(mol):
    return Descriptors.SlogP_VSA2(mol)

def getSlogP_VSA3(mol):
    return Descriptors.SlogP_VSA3(mol)

def getSlogP_VSA4(mol):
    return Descriptors.SlogP_VSA4(mol)

def getSlogP_VSA5(mol):
    return Descriptors.SlogP_VSA5(mol)

def getSlogP_VSA6(mol):
    return Descriptors.SlogP_VSA6(mol)

def getSlogP_VSA7(mol):
    return Descriptors.SlogP_VSA7(mol)

def getSlogP_VSA8(mol):
    return Descriptors.SlogP_VSA8(mol)

def getSlogP_VSA9(mol):
    return Descriptors.SlogP_VSA9(mol)

def getSlogP_VSA10(mol):
    return Descriptors.SlogP_VSA10(mol)

def getSlogP_VSA11(mol):
    return Descriptors.SlogP_VSA11(mol)

def getSlogP_VSA12(mol):
    return Descriptors.SlogP_VSA12(mol)

def getSMR_VSA1(mol):
    return Descriptors.SMR_VSA1(mol)

def getSMR_VSA2(mol):
    return Descriptors.SMR_VSA2(mol)

def getSMR_VSA3(mol):
    return Descriptors.SMR_VSA3(mol)

def getSMR_VSA4(mol):
    return Descriptors.SMR_VSA4(mol)

def getSMR_VSA5(mol):
    return Descriptors.SMR_VSA5(mol)

def getSMR_VSA6(mol):
    return Descriptors.SMR_VSA6(mol)

def getSMR_VSA7(mol):
    return Descriptors.SMR_VSA7(mol)

def getSMR_VSA8(mol):
    return Descriptors.SMR_VSA8(mol)

def getSMR_VSA9(mol):
    return Descriptors.SMR_VSA9(mol)

def getSMR_VSA10(mol):
    return Descriptors.SMR_VSA10(mol)

def getPEOE_VSA1(mol):
    return Descriptors.PEOE_VSA1(mol)

def getPEOE_VSA2(mol):
    return Descriptors.PEOE_VSA2(mol)

def getPEOE_VSA3(mol):
    return Descriptors.PEOE_VSA3(mol)

def getPEOE_VSA4(mol):
    return Descriptors.PEOE_VSA4(mol)

def getPEOE_VSA5(mol):
    return Descriptors.PEOE_VSA5(mol)

def getPEOE_VSA6(mol):
    return Descriptors.PEOE_VSA6(mol)

def getPEOE_VSA7(mol):
    return Descriptors.PEOE_VSA7(mol)

def getPEOE_VSA8(mol):
    return Descriptors.PEOE_VSA8(mol)

def getPEOE_VSA9(mol):
    return Descriptors.PEOE_VSA9(mol)

def getPEOE_VSA10(mol):
    return Descriptors.PEOE_VSA10(mol)

def getPEOE_VSA11(mol):
    return Descriptors.PEOE_VSA11(mol)

def getPEOE_VSA12(mol):
    return Descriptors.PEOE_VSA12(mol)

def getPEOE_VSA13(mol):
    return Descriptors.PEOE_VSA13(mol)

def getPEOE_VSA14(mol):
    return Descriptors.PEOE_VSA14(mol)

def getVSA_EState1(mol):
    return Descriptors.VSA_EState1(mol)

def getVSA_EState2(mol):
    return Descriptors.VSA_EState2(mol)

def getVSA_EState3(mol):
    return Descriptors.VSA_EState3(mol)

def getVSA_EState4(mol):
    return Descriptors.VSA_EState4(mol)

def getVSA_EState5(mol):
    return Descriptors.VSA_EState5(mol)

def getVSA_EState6(mol):
    return Descriptors.VSA_EState6(mol)

def getVSA_EState7(mol):
    return Descriptors.VSA_EState7(mol)

def getVSA_EState8(mol):
    return Descriptors.VSA_EState8(mol)

def getVSA_EState9(mol):
    return Descriptors.VSA_EState9(mol)

def getVSA_EState10(mol):
    return Descriptors.VSA_EState10(mol)



_moe = {"LabuteASA": getLabuteASA,
        "TPSA": getTPSA,
        "SlogP_VSA1": getSlogP_VSA1,
        "SlogP_VSA2": getSlogP_VSA2,
        "SlogP_VSA3": getSlogP_VSA3,
        "SlogP_VSA4": getSlogP_VSA4,
        "SlogP_VSA5": getSlogP_VSA5,
        "SlogP_VSA6": getSlogP_VSA6,
        "SlogP_VSA7": getSlogP_VSA7,
        "SlogP_VSA8": getSlogP_VSA8,
        "SlogP_VSA9": getSlogP_VSA9,
        "SlogP_VSA10": getSlogP_VSA10,
        "SlogP_VSA11": getSlogP_VSA11,
        "SlogP_VSA12": getSlogP_VSA12,
        "SMR_VSA1": getSMR_VSA1,
        "SMR_VSA2": getSMR_VSA2,
        "SMR_VSA3": getSMR_VSA3,
        "SMR_VSA4": getSMR_VSA4,
        "SMR_VSA5": getSMR_VSA5,
        "SMR_VSA6": getSMR_VSA6,
        "SMR_VSA7": getSMR_VSA7,
        "SMR_VSA8": getSMR_VSA8,
        "SMR_VSA9": getSMR_VSA9,
        "SMR_VSA10": getSMR_VSA10,
        "PEOE_VSA1": getPEOE_VSA1,
        "PEOE_VSA2": getPEOE_VSA2,
        "PEOE_VSA3": getPEOE_VSA3,
        "PEOE_VSA4": getPEOE_VSA4,
        "PEOE_VSA5": getPEOE_VSA5,
        "PEOE_VSA6": getPEOE_VSA6,
        "PEOE_VSA7": getPEOE_VSA7,
        "PEOE_VSA8": getPEOE_VSA8,
        "PEOE_VSA9": getPEOE_VSA9,
        "PEOE_VSA10": getPEOE_VSA10,
        "PEOE_VSA11": getPEOE_VSA11,
        "PEOE_VSA12": getPEOE_VSA12,
        "PEOE_VSA13": getPEOE_VSA13,
        "PEOE_VSA14": getPEOE_VSA14,
        "VSA_EState1": getVSA_EState1,
        "VSA_EState2": getVSA_EState2,
        "VSA_EState3": getVSA_EState3,
        "VSA_EState4": getVSA_EState4,
        "VSA_EState5": getVSA_EState5,
        "VSA_EState6": getVSA_EState6,
        "VSA_EState7": getVSA_EState7,
        "VSA_EState8": getVSA_EState8,
        "VSA_EState9": getVSA_EState9,
        "VSA_EState10": getVSA_EState10}


def GetMOE(mol):
    result={}
    for DesLabel in _moe.keys():
        result[DesLabel] = _moe[DesLabel](mol)
    return result


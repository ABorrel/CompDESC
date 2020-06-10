from rdkit.Chem.EState import EState, Fingerprinter, EState_VSA


######### fingerprint ######

def getCfragSEStatefrag(mol):
    dout = {}
    dout["Cfrag"] = {}
    dout["SEStatefrag"] = {}
    counts, sums = Fingerprinter.FingerprintMol(mol)
    for i in range(1, len(counts)+1):
        dout["Cfrag"]["Cfrag%s"%(i)] = counts[i-1]
        dout["SEStatefrag"]["SEStatefrag%s" % (i)] = sums[i - 1]
    return dout

######### based on molecule ######

def getMaxEStateIndex(mol):
    return EState.MaxEStateIndex(mol)

def getMinEStateIndex(mol):
    return EState.MinEStateIndex(mol)

def getMaxAbsEStateIndex(mol):
    return EState.MaxAbsEStateIndex(mol)

def getMinAbsEStateIndex(mol):
    return EState.MinAbsEStateIndex(mol)

######### based on the VSA #######

def getEStateVSA1(mol):
    return EState_VSA.EState_VSA1(mol)

def getEStateVSA2(mol):
    return EState_VSA.EState_VSA2(mol)

def getEStateVSA3(mol):
    return EState_VSA.EState_VSA3(mol)

def getEStateVSA4(mol):
    return EState_VSA.EState_VSA4(mol)

def getEStateVSA5(mol):
    return EState_VSA.EState_VSA5(mol)

def getEStateVSA6(mol):
    return EState_VSA.EState_VSA6(mol)

def getEStateVSA7(mol):
    return EState_VSA.EState_VSA7(mol)

def getEStateVSA8(mol):
    return EState_VSA.EState_VSA8(mol)

def getEStateVSA9(mol):
    return EState_VSA.EState_VSA9(mol)

def getEStateVSA10(mol):
    return EState_VSA.EState_VSA10(mol)

_EState = {"Cfrag1": getCfragSEStatefrag,
           "Cfrag2": getCfragSEStatefrag,
           "Cfrag3": getCfragSEStatefrag,
           "Cfrag4": getCfragSEStatefrag,
           "Cfrag5": getCfragSEStatefrag,
           "Cfrag6": getCfragSEStatefrag,
           "Cfrag7": getCfragSEStatefrag,
           "Cfrag8": getCfragSEStatefrag,
           "Cfrag9": getCfragSEStatefrag,
           "Cfrag10": getCfragSEStatefrag,
           "Cfrag11": getCfragSEStatefrag,
           "Cfrag12": getCfragSEStatefrag,
           "Cfrag13": getCfragSEStatefrag,
           "Cfrag14": getCfragSEStatefrag,
           "Cfrag15": getCfragSEStatefrag,
           "Cfrag16": getCfragSEStatefrag,
           "Cfrag17": getCfragSEStatefrag,
           "Cfrag18": getCfragSEStatefrag,
           "Cfrag19": getCfragSEStatefrag,
           "Cfrag20": getCfragSEStatefrag,
           "Cfrag21": getCfragSEStatefrag,
           "Cfrag22": getCfragSEStatefrag,
           "Cfrag23": getCfragSEStatefrag,
           "Cfrag24": getCfragSEStatefrag,
           "Cfrag25": getCfragSEStatefrag,
           "Cfrag26": getCfragSEStatefrag,
           "Cfrag27": getCfragSEStatefrag,
           "Cfrag28": getCfragSEStatefrag,
           "Cfrag29": getCfragSEStatefrag,
           "Cfrag30": getCfragSEStatefrag,
           "Cfrag31": getCfragSEStatefrag,
           "Cfrag32": getCfragSEStatefrag,
           "Cfrag33": getCfragSEStatefrag,
           "Cfrag34": getCfragSEStatefrag,
           "Cfrag35": getCfragSEStatefrag,
           "Cfrag36": getCfragSEStatefrag,
           "Cfrag37": getCfragSEStatefrag,
           "Cfrag38": getCfragSEStatefrag,
           "Cfrag39": getCfragSEStatefrag,
           "Cfrag40": getCfragSEStatefrag,
           "Cfrag41": getCfragSEStatefrag,
           "Cfrag42": getCfragSEStatefrag,
           "Cfrag43": getCfragSEStatefrag,
           "Cfrag44": getCfragSEStatefrag,
           'Cfrag45': getCfragSEStatefrag,
           "Cfrag46": getCfragSEStatefrag,
           "Cfrag47": getCfragSEStatefrag,
           "Cfrag48": getCfragSEStatefrag,
           "Cfrag49": getCfragSEStatefrag,
           "Cfrag50": getCfragSEStatefrag,
           "Cfrag51": getCfragSEStatefrag,
           "Cfrag52": getCfragSEStatefrag,
           "Cfrag53": getCfragSEStatefrag,
           "Cfrag54": getCfragSEStatefrag,
           "Cfrag55": getCfragSEStatefrag,
           "Cfrag56": getCfragSEStatefrag,
           "Cfrag57": getCfragSEStatefrag,
           "Cfrag58": getCfragSEStatefrag,
           "Cfrag59": getCfragSEStatefrag,
           "Cfrag60": getCfragSEStatefrag,
           "Cfrag61": getCfragSEStatefrag,
           "Cfrag62": getCfragSEStatefrag,
           "Cfrag63": getCfragSEStatefrag,
           "Cfrag64": getCfragSEStatefrag,
           "Cfrag65": getCfragSEStatefrag,
           "Cfrag66": getCfragSEStatefrag,
           "Cfrag67": getCfragSEStatefrag,
           "Cfrag68": getCfragSEStatefrag,
           "Cfrag69": getCfragSEStatefrag,
           "Cfrag70": getCfragSEStatefrag,
           "Cfrag71": getCfragSEStatefrag,
           "Cfrag72": getCfragSEStatefrag,
           "Cfrag73": getCfragSEStatefrag,
           "Cfrag74": getCfragSEStatefrag,
           "Cfrag75": getCfragSEStatefrag,
           "Cfrag76": getCfragSEStatefrag,
           "Cfrag77": getCfragSEStatefrag,
           "Cfrag78": getCfragSEStatefrag,
           "Cfrag79": getCfragSEStatefrag,
           "SEStatefrag1": getCfragSEStatefrag,
           "SEStatefrag2": getCfragSEStatefrag,
           "SEStatefrag3": getCfragSEStatefrag,
           "SEStatefrag4": getCfragSEStatefrag,
           "SEStatefrag5": getCfragSEStatefrag,
           "SEStatefrag6": getCfragSEStatefrag,
           "SEStatefrag7": getCfragSEStatefrag,
           "SEStatefrag8": getCfragSEStatefrag,
           "SEStatefrag9": getCfragSEStatefrag,
           "SEStatefrag10": getCfragSEStatefrag,
           "SEStatefrag11": getCfragSEStatefrag,
           "SEStatefrag12": getCfragSEStatefrag,
           "SEStatefrag13": getCfragSEStatefrag,
           "SEStatefrag14": getCfragSEStatefrag,
           "SEStatefrag15": getCfragSEStatefrag,
           "SEStatefrag16": getCfragSEStatefrag,
           "SEStatefrag17": getCfragSEStatefrag,
           "SEStatefrag18": getCfragSEStatefrag,
           "SEStatefrag19": getCfragSEStatefrag,
           "SEStatefrag20": getCfragSEStatefrag,
           "SEStatefrag21": getCfragSEStatefrag,
           "SEStatefrag22": getCfragSEStatefrag,
           "SEStatefrag23": getCfragSEStatefrag,
           "SEStatefrag24": getCfragSEStatefrag,
           "SEStatefrag25": getCfragSEStatefrag,
           "SEStatefrag26": getCfragSEStatefrag,
           "SEStatefrag27": getCfragSEStatefrag,
           "SEStatefrag28": getCfragSEStatefrag,
           "SEStatefrag29": getCfragSEStatefrag,
           "SEStatefrag30": getCfragSEStatefrag,
           "SEStatefrag31": getCfragSEStatefrag,
           "SEStatefrag32": getCfragSEStatefrag,
           "SEStatefrag33": getCfragSEStatefrag,
           "SEStatefrag34": getCfragSEStatefrag,
           "SEStatefrag35": getCfragSEStatefrag,
           "SEStatefrag36": getCfragSEStatefrag,
           "SEStatefrag37": getCfragSEStatefrag,
           "SEStatefrag38": getCfragSEStatefrag,
           "SEStatefrag39": getCfragSEStatefrag,
           "SEStatefrag40": getCfragSEStatefrag,
           "SEStatefrag41": getCfragSEStatefrag,
           "SEStatefrag42": getCfragSEStatefrag,
           "SEStatefrag43": getCfragSEStatefrag,
           "SEStatefrag44": getCfragSEStatefrag,
           "SEStatefrag45": getCfragSEStatefrag,
           "SEStatefrag46": getCfragSEStatefrag,
           "SEStatefrag47": getCfragSEStatefrag,
           "SEStatefrag48": getCfragSEStatefrag,
           "SEStatefrag49": getCfragSEStatefrag,
           "SEStatefrag50": getCfragSEStatefrag,
           "SEStatefrag51": getCfragSEStatefrag,
           "SEStatefrag52": getCfragSEStatefrag,
           "SEStatefrag53": getCfragSEStatefrag,
           "SEStatefrag54": getCfragSEStatefrag,
           "SEStatefrag55": getCfragSEStatefrag,
           "SEStatefrag56": getCfragSEStatefrag,
           "SEStatefrag57": getCfragSEStatefrag,
           "SEStatefrag58": getCfragSEStatefrag,
           "SEStatefrag59": getCfragSEStatefrag,
           "SEStatefrag60": getCfragSEStatefrag,
           "SEStatefrag61": getCfragSEStatefrag,
           "SEStatefrag62": getCfragSEStatefrag,
           "SEStatefrag63": getCfragSEStatefrag,
           "SEStatefrag64": getCfragSEStatefrag,
           "SEStatefrag65": getCfragSEStatefrag,
           "SEStatefrag66": getCfragSEStatefrag,
           "SEStatefrag67": getCfragSEStatefrag,
           "SEStatefrag68": getCfragSEStatefrag,
           "SEStatefrag69": getCfragSEStatefrag,
           "SEStatefrag70": getCfragSEStatefrag,
           "SEStatefrag71": getCfragSEStatefrag,
           "SEStatefrag72": getCfragSEStatefrag,
           "SEStatefrag73": getCfragSEStatefrag,
           "SEStatefrag74": getCfragSEStatefrag,
           "SEStatefrag75": getCfragSEStatefrag,
           "SEStatefrag76": getCfragSEStatefrag,
           "SEStatefrag77": getCfragSEStatefrag,
           "SEStatefrag78": getCfragSEStatefrag,
           "SEStatefrag79": getCfragSEStatefrag,
           "MaxEStateIndex": getMaxEStateIndex,
           "MinEStateIndex": getMinEStateIndex,
           "MaxAbsEStateIndex": getMaxAbsEStateIndex,
           "MinAbsEStateIndex": getMinAbsEStateIndex,
           "EStateVSA1": getEStateVSA1,
           "EStateVSA2": getEStateVSA2,
           "EStateVSA3": getEStateVSA3,
           "EStateVSA4": getEStateVSA4,
           "EStateVSA5": getEStateVSA5,
           "EStateVSA6": getEStateVSA6,
           "EStateVSA7": getEStateVSA7,
           "EStateVSA8": getEStateVSA8,
           "EStateVSA9": getEStateVSA9,
           "EStateVSA10": getEStateVSA10}


def GetEState(mol):
    dresult = {}
    # fingerprint
    dfingerEState = getCfragSEStatefrag(mol)
    dresult.update(dfingerEState["Cfrag"])
    dresult.update(dfingerEState["SEStatefrag"])
    # by atom
    dresult["MaxEStateIndex"] = getMaxEStateIndex(mol)
    dresult["MinEStateIndex"] = getMinEStateIndex(mol)
    dresult["MaxAbsEStateIndex"] = getMaxAbsEStateIndex(mol)
    dresult["MinAbsEStateIndex"] = getMinAbsEStateIndex(mol)
    # VSA
    dresult["EStateVSA1"] = getEStateVSA1(mol)
    dresult["EStateVSA2"] = getEStateVSA2(mol)
    dresult["EStateVSA3"] = getEStateVSA3(mol)
    dresult["EStateVSA4"] = getEStateVSA4(mol)
    dresult["EStateVSA5"] = getEStateVSA5(mol)
    dresult["EStateVSA6"] = getEStateVSA6(mol)
    dresult["EStateVSA7"] = getEStateVSA7(mol)
    dresult["EStateVSA8"] = getEStateVSA8(mol)
    dresult["EStateVSA9"] = getEStateVSA9(mol)
    dresult["EStateVSA10"] = getEStateVSA10(mol)
    return dresult

#from rdkit import Chem
#mol = Chem.MolFromSmiles("OC(=NC1=NC2=C(N=CN2C2CC(OC(=O)C3=CC=CC=C3)C(COC(=O)C3=CC=CC=C3)O2)C(O)=N1)C1=CC=CC=C1")
#d = GetEState(mol)
#print(d)



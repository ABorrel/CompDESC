from rdkit.Chem import Descriptors


def getpkappa1(mol):
    P1 = float(mol.GetNumBonds(1))
    A = float(mol.GetNumHeavyAtoms())
    num = A * ((A-1) ** 2)
    denom = P1 **2
    if denom != 0.0 :
        kappa = num / denom
    else:
        kappa = 0.0
    return kappa

def getpkappa2(mol):
    P2 = float(mol.GetNumBonds(2))
    A = float(mol.GetNumHeavyAtoms())
    num = (A - 1) * ((A-2) ** 2)
    denom = P2 **2
    if denom != 0.0 :
        kappa = num / denom
    else:
        kappa = 0.0
    return kappa

def getpkappa3(mol):
    P3 = float(mol.GetNumBonds(3))
    A = float(mol.GetNumHeavyAtoms())
    num = (A - 1) * ((A-3) ** 2)
    denom = P3 **2
    if denom != 0.0 :
        kappa = num / denom
    else:
        kappa = 0.0
    return kappa

def getskappa1(mol):
    return Descriptors.Kappa1(mol)

def getskappa2(mol):
    return Descriptors.Kappa2(mol)

def getskappa3(mol):
    return Descriptors.Kappa3(mol)

def getphi(mol):
    """
    #################################################################
    Calculation of Kier molecular flexibility index
    ---->phi
    #################################################################
    """
    skappa1 = getskappa1(mol)
    skappa2 = getskappa2(mol)
    A = float(mol.GetNumHeavyAtoms())
    phi = (skappa1*skappa2)/A
    return phi

def getHallKierAlpha(mol):
    return Descriptors.HallKierAlpha(mol)


_kappa = {'pkappa1': getpkappa1,
          'pkappa2': getpkappa2,
          'pkappa3': getpkappa3,
          'skappa1': getskappa1,
          'skappa2': getskappa2,
          'skappa3': getskappa3,
          'phi': getphi,
          'HallKierAlpha': getHallKierAlpha}

def GetKappa(mol):
    dresult = {}
    for DesLabel in _kappa.keys():
        dresult[DesLabel] = round(_kappa[DesLabel](mol), 6)
    return dresult

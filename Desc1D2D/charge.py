# Should put in a class to not repeat the charge position
from rdkit import Chem
from rdkit.Chem import rdPartialCharges as GMCharge, Descriptors

import numpy
import math
iter_step=12


def getElementAtomCharges(mol, AtomicNum):
    """
    #################################################################
    Most positive charge on atom with atomic number equal to n
    #################################################################
    """
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    lcharge = []
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNum:
            charge = float(atom.GetProp('_GasteigerCharge'))
            if not math.isnan(charge) and not charge == numpy.inf:
                lcharge.append(float(atom.GetProp('_GasteigerCharge')))
        if AtomicNum == "All":
            charge = float(atom.GetProp('_GasteigerCharge'))
            if not math.isnan(charge) and not charge == numpy.inf:
                lcharge.append(float(atom.GetProp('_GasteigerCharge')))
    if lcharge == []:
        lcharge = [0.0]
    return lcharge

def getElementAtomSumSquareCharge(mol, AtomicNum):
    """
    #################################################################
    Ths sum of square Charges on all atoms with atomicnumber equal to n
    #################################################################
    """
    lcharges = getElementAtomCharges(mol, AtomicNum)
    if lcharges == []:
        return 0.0
    else:
        return round(sum(numpy.square(lcharges)), 6)


#################################################

def getQHmax(mol):
    return round(max(getElementAtomCharges(mol, AtomicNum=1)), 6)

def getQCmax(mol):
    return round(max(getElementAtomCharges(mol, AtomicNum=6)), 6)

def getQNmax(mol):
    return round(max(getElementAtomCharges(mol, AtomicNum=7)), 6)

def getQOmax(mol):
    return round(max(getElementAtomCharges(mol, AtomicNum=8)), 6)

def getQHmin(mol):
    return round(min(getElementAtomCharges(mol, AtomicNum=1)), 6)

def getQCmin(mol):
    return round(min(getElementAtomCharges(mol, AtomicNum=6)), 6)

def getQNmin(mol):
    return round(min(getElementAtomCharges(mol, AtomicNum=7)), 6)

def getQOmin(mol):
    return round(min(getElementAtomCharges(mol, AtomicNum=8)), 6)

def getQmin(mol):
    return round(min(getElementAtomCharges(mol, AtomicNum="All")), 6)

def getQmax(mol):
    return round(max(getElementAtomCharges(mol, AtomicNum="All")), 6)

def getQss(mol):
    return getElementAtomSumSquareCharge(mol, "All")

def getQCss(mol):
    return getElementAtomSumSquareCharge(mol, 6)

def getQNss(mol):
    return getElementAtomSumSquareCharge(mol, 7)

def getQOss(mol):
    return getElementAtomSumSquareCharge(mol, 8)

def getQHss(mol):
    return getElementAtomSumSquareCharge(mol, 1)

def getTpc(mol):
    """
    #################################################################
    The total postive charge
    #################################################################
    """
    lcharges = getElementAtomCharges(mol, "All")
    lp = []
    for charge in lcharges:
        if charge > 0.0:
            lp.append(charge)
    if lp == []:
        return 0.0
    else:
        return round(sum(lp),6)

def getMpc(mol):
    """
    #################################################################
    The average postive charge
    #################################################################
    """
    lcharges = getElementAtomCharges(mol, "All")
    lp = []
    for charge in lcharges:
        if charge > 0.0:
            lp.append(charge)
    if lp == []:
        return 0.0
    else:
        return round(numpy.mean(lp), 6)

def getTnc(mol):
    """
    #################################################################
    The total negative charge
    #################################################################
    """
    lcharges = getElementAtomCharges(mol, "All")
    ln = []
    for charge in lcharges:
        if charge < 0.0:
            ln.append(charge)
    if ln == []:
        return 0.0
    else:
        return round(sum(ln),6)

def getMnc(mol):
    """
    #################################################################
    The average postive charge
    #################################################################
    """
    lcharges = getElementAtomCharges(mol, "All")
    ln = []
    for charge in lcharges:
        if charge < 0.0:
            ln.append(charge)
    if ln == []:
        return 0.0
    else:
        return round(numpy.mean(ln), 6)

def getTac(mol):
    lcharges = getElementAtomCharges(mol, "All")
    lcharges = numpy.absolute(lcharges)
    if len(lcharges) == 0:
        return 0.0
    else:
        return round(sum(lcharges), 6)

def getMac(mol):
    lcharges = getElementAtomCharges(mol, "All")
    lcharges = numpy.absolute(lcharges)
    if len(lcharges) == 0:
        return 0.0
    else:
        return round(numpy.mean(lcharges), 6)

def getRpc(mol):
    lcharges = getElementAtomCharges(mol, "All")
    lp = []

    for charge in lcharges:
        if charge > 0.0:
            lp.append(charge)
    if lp == []:
        return 0.0
    else:
        return round(max(lp)/sum(lp), 6)

def getRnc(mol):
    lcharges = getElementAtomCharges(mol, "All")
    ln = []

    for charge in lcharges:
        if charge < 0.0:
            ln.append(charge)
    if ln == []:
        return 0.0
    else:
        return round(min(ln)/sum(ln), 6)

def getLDI(mol):
    """
    #################################################################
    Calculation of local dipole index (D)
    #################################################################
    """
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
        charge = float(atom.GetProp('_GasteigerCharge'))
        if not math.isnan(charge) and not charge == numpy.inf:
            res.append(charge)
        else:
            res.append(0.0)
    cc = [numpy.absolute(res[x.GetBeginAtom().GetIdx()]-res[x.GetEndAtom().GetIdx()]) for x in Hmol.GetBonds()]
    B = len(Hmol.GetBonds())
    if B == 0:
        return 0.0
    return round(sum(cc)/B, 6)
        
def getSPP(mol):
    """
    #################################################################
    Calculation of submolecular polarity parameter(SPP)
    #################################################################
    """
    return round(getQmax(mol)-getQmin(mol), 6)

def getNumValenceElectrons(mol):
    Hmol = Chem.AddHs(mol)
    return Descriptors.NumValenceElectrons(Hmol)

def getNumRadicalElectrons(mol):
    Hmol = Chem.AddHs(mol)
    return Descriptors.NumRadicalElectrons(Hmol)

def getMaxAbsPartialCharge(mol):
    Hmol = Chem.AddHs(mol)
    abscharge = Descriptors.MaxAbsPartialCharge(Hmol)
    if abscharge == numpy.inf or math.isnan(abscharge) == True:
        return 0.0
    else: 
        return abscharge

def getMinAbsPartialCharge(mol):
    Hmol = Chem.AddHs(mol)
    absmin = Descriptors.MinAbsPartialCharge(Hmol)
    if absmin == numpy.inf or math.isnan(absmin) == True :
        return 0.0
    else:
        return absmin


_charge={'SPP':getSPP,
        'LDI':getLDI,
        'Rnc':getRnc,
        'Rpc':getRpc,
        'Mac':getMac,
        'Tac':getTac,
        'Mnc':getMnc,
        'Tnc':getTnc,
        'Mpc':getMpc,
        'Tpc':getTpc,
        'Qss':getQss,
        'QOss':getQOss,
        'QNss':getQNss,
        'QCss':getQCss,
        'QHss':getQHss,
        'Qmin':getQmin,
        'Qmax':getQmax,
        'QOmin':getQOmin,
        'QNmin':getQNmin,
        'QCmin':getQCmin,
        'QHmin':getQHmin,
        'QOmax':getQOmax,
        'QNmax':getQNmax,
        'QCmax':getQCmax,
        'QHmax':getQHmax,
        'NumValenceElectrons':getNumValenceElectrons,
        'NumRadicalElectrons':getNumRadicalElectrons,
        'MaxAbsPartialCharge':getMaxAbsPartialCharge,
        'MinAbsPartialCharge':getMinAbsPartialCharge}


def GetCharge(mol):
    result={}
    for DesLabel in _charge.keys():
        result[DesLabel] =_charge[DesLabel](mol)
    return result


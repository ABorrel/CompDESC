from rdkit import Chem
from rdkit.Chem import rdchem, Descriptors
import numpy


periodicTable = rdchem.GetPeriodicTable()

def getChinp(mol,NumPath=2):
    """
    #################################################################
    Calculation of molecular connectivity chi index for path order n
    #################################################################
    """
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    for path in Chem.FindAllPathsOfLengthN(mol,NumPath+1,useBonds=0):
        cAccum=1.0
        for idx in path:
            cAccum *= deltas[idx]
        if cAccum:
            accum += 1./numpy.sqrt(cAccum)
    return accum


def getChinch(mol, NumCycle=3):
    """
    #################################################################
    Calculation of molecular connectivity chi index for cycles of n
    #################################################################
    """
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    for tup in mol.GetRingInfo().AtomRings():
        cAccum = 1.0
        if len(tup) == NumCycle:
            for idx in tup:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1. / numpy.sqrt(cAccum)
    return accum

def getHKDeltas(mol, skipHs=1):
    """
    #################################################################
    Calculation of modified delta value for a molecule
    #################################################################
    """
    global periodicTable
    res = []
    for atom in mol.GetAtoms():
        n = atom.GetAtomicNum()
        if n > 1:
            nV = periodicTable.GetNOuterElecs(n)
            nHs = atom.GetTotalNumHs()
            if n < 10:
                res.append(float(nV - nHs))
            else:
                res.append(float(nV - nHs) / float(n - nV - 1))
        elif not skipHs:
            res.append(0.0)
    return res

def getAtomHKDeltas(atom,skipHs=0):
    """
    #################################################################
    *Internal Use Only*
    Calculation of modified delta value for a molecule
    #################################################################
    """
    global periodicTable
    res=[]
    n=atom.GetAtomicNum()
    if n > 1:
        nV=periodicTable.GetNOuterElecs(n)
        nHs=atom.GetTotalNumHs()
        if n<10:
            res.append(float(nV-nHs))
        else:
            res.append(float(nV-nHs)/float(n-nV-1))
    elif not skipHs:
        res.append(0.0)
    return res

def getChivnp(mol, NumPath=1):
    """#################################################################
    Calculation of valence molecular connectivity chi index for path order 1
    #################################################################
    """
    accum = 0.0
    deltas = getHKDeltas(mol, skipHs=0)
    for path in Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0):
        cAccum = 1.0
        for idx in path:
            cAccum *= deltas[idx]
        if cAccum:
            accum += 1. / numpy.sqrt(cAccum)
    return accum

def getChivnch(mol, NumCyc=3):
    """
    #################################################################
    Calculation of valence molecular connectivity chi index for cycles of n
    #################################################################
    """
    accum=0.0
    deltas=getHKDeltas(mol,skipHs=0)
    for tup in mol.GetRingInfo().AtomRings():
        cAccum=1.0
        if len(tup)==NumCyc:
            for idx in tup:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1./numpy.sqrt(cAccum)
    return accum
################################################################

def getChi0(mol):
    return Descriptors.Chi0(mol)

def getChi1(mol):
    return Descriptors.Chi1(mol)

def getmChi1(mol):
    """
    #################################################################
    Calculation of mean chi1 (Randic) connectivity index.
    ---->mchi1
    #################################################################
    """
    cc = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc,'d')
    res = numpy.mean(numpy.sqrt(1./cc))
    return res

def getChi2(mol):
    return getChinp(mol,NumPath=2)

def getChi3(mol):
    return getChinp(mol,NumPath=3)

def getChi4(mol):
    return getChinp(mol,NumPath=4)

def getChi5(mol):
    return getChinp(mol,NumPath=5)

def getChi6(mol):
    return getChinp(mol,NumPath=6)

def getChi7(mol):
    return getChinp(mol,NumPath=7)

def getChi8(mol):
    return getChinp(mol,NumPath=8)

def getChi9(mol):
    return getChinp(mol,NumPath=9)

def getChi10(mol):
    return getChinp(mol,NumPath=10)


def getChi3c(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas=[mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum

def getChi4c(mol):
    accum=0.0
    patt=Chem.MolFromSmarts('*~*(~*)(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas=[mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,float)
            accum=accum + 1./numpy.sqrt(deltas1.prod())
    return accum
    

def getChi4pc(mol):
    accum=0.0
    patt=Chem.MolFromSmarts('*~*(~*)~*~*')
    HPatt=mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas=[mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum

def getChi3ch(mol):
    return getChinch(mol,NumCycle=3)

def getChi4ch(mol):
    return getChinch(mol,NumCycle=4)

def getChi5ch(mol):
    return getChinch(mol,NumCycle=5)

def getChi6ch(mol):
    return getChinch(mol,NumCycle=6)

def getChiv0(mol):
    deltas=getHKDeltas(mol,skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    res=sum(numpy.sqrt(1./deltas))
    return res

def getChiv1(mol):
    return getChivnp(mol,NumPath=1)
    
def getChiv2(mol):
    return getChivnp(mol,NumPath=2)

def getChiv3(mol):
    return getChivnp(mol, NumPath=3)

def getChiv4(mol):
    return getChivnp(mol, NumPath=4)

def getChiv5(mol):
    return getChivnp(mol, NumPath=5)

def getChiv6(mol):
    return getChivnp(mol, NumPath=6)

def getChiv7(mol):
    return getChivnp(mol, NumPath=7)

def getChiv8(mol):
    return getChivnp(mol, NumPath=8)

def getChiv9(mol):
    return getChivnp(mol, NumPath=9)

def getChiv10(mol):
    return getChivnp(mol, NumPath=10)

def getdchi0(mol):
    """
    #################################################################
    Calculation of the difference between chi0v and chi0
    #################################################################
    """
    return abs(getChiv0(mol) - getChi0(mol))
   
    
def getdchi1(mol):
    return abs(getChiv1(mol)-getChi1(mol))

def getdchi2(mol):
    return abs(getChiv2(mol)-getChi2(mol))

def getdchi3(mol):
    return abs(getChiv3(mol)-getChi3(mol))

def getdchi4(mol):
    return abs(getChiv4(mol)-getChi4(mol))

def getChiv3c(mol):
    accum=0.0
    patt=Chem.MolFromSmarts('*~*(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas=[getAtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,float)
            den = numpy.sqrt(deltas1.prod())
            if den != 0.0:
                accum = accum + 1./den
    return accum

def getChiv4c(mol):
    accum=0.0
    patt=Chem.MolFromSmarts('*~*(~*)(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[getAtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,float)
            den = numpy.sqrt(deltas1.prod())
            if den != 0.0:
                accum = accum+1./den
    return accum

def getChiv4pc(mol):
    """
    #################################################################
    Calculation of valence molecular connectivity chi index for 
    #################################################################
    """
    accum=0.0
    patt=Chem.MolFromSmarts('*~*(~*)~*~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[getAtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,float)
            den = numpy.sqrt(deltas1.prod())
            if den != 0.0:
                accum = accum+1./den
    return accum

def getChiv3ch(mol):
    """
    Chiv3ch related to ring  3
    """
    return getChivnch(mol, 3)

def getChiv4ch(mol):
    """
    Chiv4ch related to ring  4
    """
    return getChivnch(mol, 4)

def getChiv5ch(mol):
    """
    Chiv5ch related to ring  5
    """
    return getChivnch(mol, 5)

def getChiv6ch(mol):
    """
    Chiv6h related to ring  6
    """
    return getChivnch(mol, 6)


def getknotp(mol):
    """
    #################################################################
    Calculation of the difference between chi3c and chi4pc
    #################################################################
    """
    return abs(getChi3c(mol)-getChi4pc(mol))


def getknotpv(mol):
    """
    #################################################################
    Calculation of the difference between chiv3c and chiv4pc
    ---->knotpv
    #################################################################
    """
    chiv3 = getChiv3c(mol)
    chiv4pc = getChiv4pc(mol)
    return abs(getChiv3c(mol) - getChiv4pc(mol))

def getChi0n(mol):
    return Chem.GraphDescriptors.Chi0n(mol)

def getChi1n(mol):
    return Chem.GraphDescriptors.Chi1n(mol)

def getChi2n(mol):
    return Chem.GraphDescriptors.Chi2n(mol)

def getChi3n(mol):
    return Chem.GraphDescriptors.Chi3n(mol)

def getChi4n(mol):
    return Chem.GraphDescriptors.Chi4n(mol)



_connectivity={"Chi0": getChi0,
               "Chi1": getChi1,
               "mChi1": getmChi1,
               "Chi2": getChi2,
               "Chi3": getChi3,
               "Chi4": getChi4,
               "Chi5": getChi5,
               "Chi6": getChi6,
               "Chi7": getChi7,
               "Chi8": getChi8,
               "Chi9": getChi9,
               "Chi10": getChi10,
               "Chi3c": getChi3c,
               "Chi4c": getChi4c,
               "Chi4pc": getChi4pc,
               "Chi3ch": getChi3ch,
               "Chi4ch": getChi4ch,
               "Chi5ch": getChi5ch,
               "Chi6ch": getChi6ch,
               "Chiv0": getChiv0,
               "Chiv1": getChiv1,
               "Chiv2": getChiv2,
               "Chiv3": getChiv3,
               "Chiv4": getChiv4,
               "Chiv5": getChiv5,
               "Chiv6": getChiv6,
               "Chiv7": getChiv7,
               "Chiv8": getChiv8,
               "Chiv9": getChiv9,
               "Chiv10": getChiv10,
               "dchi0": getdchi0,
               "dchi1": getdchi1,
               "dchi2": getdchi2,
               "dchi3": getdchi3,
               "dchi4": getdchi4,
               "Chiv3c": getChiv3c,
               "Chiv4c": getChiv4c,
               "Chiv4pc": getChiv4pc,
               "Chiv3ch": getChiv3ch,
               "Chiv4ch": getChiv4ch,
               "Chiv5ch": getChiv5ch,
               "Chiv6ch": getChiv6ch,
               "knotp": getknotp,
               "knotpv": getknotpv,
               "Chi1n": getChi1n,
               "Chi2n": getChi2n,
               "Chi3n": getChi3n,
               "Chi4n": getChi4n}



def GetConnectivity(mol):
    """
    #################################################################
    Get the dictionary of connectivity descriptors for given moelcule mol
    #################################################################
    """
    dresult={}
    for DesLabel in _connectivity.keys():
        dresult[DesLabel]=round(_connectivity[DesLabel](mol),3)
    return dresult
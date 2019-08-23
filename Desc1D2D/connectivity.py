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


################################################################

def getChi0(mol):
    return Descriptors.Chi0(mol)

def getChi1(mol):
    return Descriptors.Chi0(mol)

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
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum

def getChi4c(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas=[mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum
    

def getChi4pc(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)~*~*')
    HPatt=mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas=[mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
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













def CalculateDeltaChi0(mol):
    """
    #################################################################
    Calculation of the difference between chi0v and chi0
    
    ---->dchi0
    
    Usage:
        
        result=CalculateDeltaChi0(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return abs(CalculateChiv0(mol)-CalculateChi0(mol))
   
    
def CalculateDeltaChi1(mol):
    """
    #################################################################
    Calculation of the difference between chi1v and chi1
    
    ---->dchi1
    
    Usage:
        
        result=CalculateDeltaChi1(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return abs(CalculateChiv1(mol)-CalculateChi1(mol))


def CalculateDeltaChi2(mol):
    """
    #################################################################
    Calculation of the difference between chi2v and chi2
    
    ---->dchi2
    
    Usage:
        
        result=CalculateDeltaChi2(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return abs(_CalculateChivnp(mol,NumPath=2)-_CalculateChinp(mol,NumPath=2))


def CalculateDeltaChi3(mol):
    """
    #################################################################
    Calculation of the difference between chi3v and chi3
    
    ---->dchi3

    Usage:
        
        result=CalculateDeltaChi3(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return abs(_CalculateChivnp(mol,NumPath=3)-_CalculateChinp(mol,NumPath=3))


def CalculateDeltaChi4(mol):
    """
    #################################################################
    Calculation of the difference between chi4v and chi4
    
    ---->dchi4

    Usage:
        
        result=CalculateDeltaChi4(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return abs(_CalculateChivnp(mol,NumPath=4)-_CalculateChinp(mol,NumPath=4))


def _AtomHKDeltas(atom,skipHs=0):
    """
    #################################################################
    *Internal Use Only*
    
    Calculation of modified delta value for a molecule
    #################################################################
    """
    global periodicTable
    res=[]
    n=atom.GetAtomicNum()
    if n>1:
        nV=periodicTable.GetNOuterElecs(n)
        nHs=atom.GetTotalNumHs()
        if n<10:
            res.append(float(nV-nHs))
        else:
            res.append(float(nV-nHs)/float(n-nV-1))
    elif not skipHs:
        res.append(0.0)
    return res


def CalculateChiv3c(mol):
    """
    #################################################################
    Calculation of valence molecular connectivity chi index for cluster
    
    ---->Chiv3c

    Usage:
        
        result=CalculateChiv3c(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[_AtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum

def CalculateChiv4c(mol):
    """
    #################################################################
    Calculation of valence molecular connectivity chi index for cluster
    
    ---->Chiv4c

    Usage:
        
        result=CalculateChiv4c(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[_AtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum

def CalculateChiv4pc(mol):
    """
    #################################################################
    Calculation of valence molecular connectivity chi index for 
    
    path/cluster
    
    ---->Chiv4pc
    
    Usage:
        
        result=CalculateChiv4pc(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)~*~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[_AtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum

def CalculateDeltaChiv3c4pc(mol):
    """
    #################################################################
    Calculation of the difference between chiv3c and chiv4pc
    
    ---->knotpv

    Usage:
        
        result=CalculateDeltaChiv3c4pc(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return abs(CalculateChiv3c(mol)-CalculateChiv4pc(mol))

def _CalculateChivnch(mol,NumCyc=3):
    """
    #################################################################
    **Internal used only**
    
    Calculation of valence molecular connectivity chi index for cycles of n
    #################################################################
    """
    accum=0.0
    deltas=_HKDeltas(mol,skipHs=0)
    for tup in mol.GetRingInfo().AtomRings():
        cAccum=1.0
        if len(tup)==NumCyc:
            for idx in tup:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1./numpy.sqrt(cAccum)

    return accum

def CalculateChiv3ch(mol):
    """
    #################################################################
    Calculation of valence molecular connectivity chi index 
    
    for cycles of 3
    
    ---->Chiv3ch

    Usage:
        
        result=CalculateChiv3ch(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return _CalculateChivnch(mol,NumCyc=3)



def CalculateChiv4ch(mol):
    """
    #################################################################
    Calculation of valence molecular connectivity chi index for 
    
    cycles of 4
    
    ---->Chiv4ch

    Usage:
        
        result=CalculateChiv4ch(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return _CalculateChivnch(mol,NumCyc=4)


def CalculateChiv5ch(mol):
    """
    #################################################################
    Calculation of valence molecular connectivity chi index for 
    
    cycles of 5
    
    ---->Chiv5ch

    Usage:
        
        result=CalculateChiv5ch(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return _CalculateChivnch(mol,NumCyc=5)



def CalculateChiv6ch(mol):
    
    """
    #################################################################
    Calculation of valence molecular connectivity chi index for
    
    cycles of 6
    
    ---->Chiv6ch

    Usage:
        
        result=CalculateChiv6ch(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value
    #################################################################
    """
    return _CalculateChivnch(mol,NumCyc=6)

def getknotp(mol):
    """
    #################################################################
    Calculation of the difference between chi3c and chi4pc
    #################################################################
    """
    return abs(getChi3c(mol)-getChi4pc(mol))





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
               "Chi0n": getChi0n,
               "Chi1n": getChi1n,
               "Chi2n": getChi2n,
               "Chi3n": getChi3n,
               "Chi4n": getChi4n}

    'Chi0':CalculateChi0,
                       'Chi1':CalculateChi1,
                       'mChi1':CalculateMeanRandic,
                       'Chi2':CalculateChi2,
                       'Chi3':CalculateChi3p,
                       'Chi4':CalculateChi4p,
                       'Chi5':CalculateChi5p,
                       'Chi6':CalculateChi6p,
                       'Chi7':CalculateChi7p,
                       'Chi8':CalculateChi8p,
                       'Chi9':CalculateChi9p,
                       'Chi10':CalculateChi10p,
                       'Chi3c':CalculateChi3c,
                       'Chi4c':CalculateChi4c,
                       'Chi4pc':CalculateChi4pc,
                       'Chi3ch':CalculateChi3ch,
                       'Chi4ch':CalculateChi4ch,
                       'Chi5ch':CalculateChi5ch,
                       'Chi6ch':CalculateChi6ch,
                       'knotp':CalculateDeltaChi3c4pc,
                       'Chiv0':CalculateChiv0,
                      'Chiv1':CalculateChiv1,
                      'Chiv2':CalculateChiv2,
                      'Chiv3':CalculateChiv3p,
                      'Chiv4':CalculateChiv4p,
                       'Chiv5':CalculateChiv5p,
                       'Chiv6':CalculateChiv6p,
                       'Chiv7':CalculateChiv7p,
                       'Chiv8':CalculateChiv8p,
                       'Chiv9':CalculateChiv9p,
                       'Chiv10':CalculateChiv10p,
                       'dchi0':CalculateDeltaChi0,
                       'dchi1':CalculateDeltaChi1,
                       'dchi2':CalculateDeltaChi2,
                       'dchi3':CalculateDeltaChi3,
                       'dchi4':CalculateDeltaChi4,
                       'Chiv3c':CalculateChiv3c,
                       'Chiv4c':CalculateChiv4c,
                       'Chiv4pc':CalculateChiv4pc,
                       'Chiv3ch':CalculateChiv3ch,
                       'Chiv4ch':CalculateChiv4ch,
                       'Chiv5ch':CalculateChiv5ch,
                       'Chiv6ch':CalculateChiv6ch,
                       'knotpv':CalculateDeltaChiv3c4pc
    }




def GetConnectivity(mol):
    """
    #################################################################
    Get the dictionary of connectivity descriptors for given moelcule mol
    
    Usage:
        
        result=GetConnectivity(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a dict form containing all connectivity indices
    #################################################################
    """
    result={}
    for DesLabel in _connectivity.keys():
        result[DesLabel]=round(_connectivity[DesLabel](mol),3)
    return result


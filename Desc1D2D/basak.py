from rdkit import Chem
import numpy
import copy


def getEntropy(Probability):
    """
    #################################################################
    Calculation of entropy (Information content) for probability given
    #################################################################
    """
    res = 0.0
    for i in Probability:
        if i != 0:
            res = res - i * numpy.log2(i)

    return res


def getICn(mol, NumPath=1):
    """
    #################################################################
    Obtain the information content with order n proposed by Basak
    #################################################################
    """
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    TotalPath = Chem.FindAllPathsOfLengthN(Hmol, NumPath, useBonds=0, useHs=1)
    if len(TotalPath) == 0:
        BasakIC = 0.0
    else:
        IC = {}
        for i in range(nAtoms):
            temp = []
            at = Hmol.GetAtomWithIdx(i)
            temp.append(at.GetAtomicNum())
            for index in TotalPath:
                if i == index[0]:
                    temp = temp + [Hmol.GetAtomWithIdx(kk).GetAtomicNum() for kk in index[1:]]
                if i == index[-1]:
                    cds = list(index)
                    cds.reverse()
                    temp = temp + [Hmol.GetAtomWithIdx(kk).GetAtomicNum() for kk in cds[1:]]

            IC[str(i)] = temp
        cds = []
        for value in IC.values():
            value.sort()
            cds.append(value)
        kkk = list(range(len(cds)))
        aaa = copy.deepcopy(kkk)
        res = []
        for i in aaa:
            if i in kkk:
                jishu = 0.0
                kong = []
                temp1 = cds[i]
                for j in aaa:
                    if cds[j] == temp1:
                        jishu = jishu + 1
                        kong.append(j)
                for ks in kong:
                    kkk.remove(ks)
                res.append(jishu)

        BasakIC = getEntropy(numpy.array(res, float) / sum(res))

    return BasakIC


def getSICn(mol, n):
    """
    #################################################################
    Obtain the structural information content with order n
    proposed by Basak.
    #################################################################
    """
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = getICn(mol, n)
    if nAtoms <= 1:
        BasakSIC = 0.0
    else:
        BasakSIC = IC / numpy.log2(nAtoms)
    return BasakSIC


def getCICn(mol, n):
    """
    #################################################################
    Obtain the complementary information content with order 1 proposed
    by Basak.
    #################################################################
    """
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = getICn(mol, n)
    if nAtoms <= 1:
        BasakCIC = 0.0
    else:
        BasakCIC = numpy.log2(nAtoms) - IC
    return BasakCIC

############################################################################

def getIC0(mol):
    """
    #################################################################
    Obtain the information content with order 0 proposed by Basak
    #################################################################
    """

    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC=[]
    for i in range(nAtoms):
        at = Hmol.GetAtomWithIdx(i)
        IC.append(at.GetAtomicNum())
    Unique=numpy.unique(IC)
    NAtomType=len(Unique)
    NTAtomType=numpy.zeros(NAtomType,float)
    for i in range(NAtomType):
        NTAtomType[i]=IC.count(Unique[i])

    if nAtoms != 0:
        #print sum(NTAtomType/nAtoms)
        BasakIC=getEntropy(NTAtomType/nAtoms)
    else:
        BasakIC = 0.0
        
    return BasakIC
        

def getSIC0(mol):
    """
    #################################################################
    Obtain the structural information content with order 0 
    #################################################################
    """

    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = getIC0(mol)
    if nAtoms <= 1:
        BasakSIC = 0.0
    else:
        BasakSIC = IC/numpy.log2(nAtoms)
    return BasakSIC
    

def getCIC0(mol):
    """
    #################################################################
    Obtain the complementary information content with order 0 
    #################################################################
    """
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = getIC0(mol)
    if nAtoms<=1:
        BasakCIC=0.0
    else:
        BasakCIC=numpy.log2(nAtoms)-IC
    return BasakCIC

def getIC1(mol):
    return getICn(mol, NumPath=2)

def getIC2(mol):
    return getICn(mol, NumPath=3)

def getIC3(mol):
    return getICn(mol, NumPath=4)

def getIC4(mol):
    return getICn(mol, NumPath=5)

def getIC5(mol):
    return getICn(mol, NumPath=6)

def getIC6(mol):
    return getICn(mol, NumPath=7)

def getSIC1(mol):
    return getSICn(mol, 2)

def getSIC2(mol):
    return getSICn(mol, 3)

def getSIC3(mol):
    return getSICn(mol, 4)

def getSIC4(mol):
    return getSICn(mol, 5)

def getSIC5(mol):
    return getSICn(mol, 6)

def getSIC6(mol):
    return getSICn(mol, 7)

def getCIC1(mol):
    return getCICn(mol, 2)

def getCIC2(mol):
    return getCICn(mol, 3)

def getCIC3(mol):
    return getCICn(mol, 4)

def getCIC4(mol):
    return getCICn(mol, 5)

def getCIC5(mol):
    return getCICn(mol, 6)

def getCIC6(mol):
    return getCICn(mol, 7)

_basak={'CIC0':getCIC0,
        'CIC1':getCIC1,
        'CIC2':getCIC2,
        'CIC3':getCIC3,
        'CIC4':getCIC4,
        'CIC5':getCIC5,
        'CIC6':getCIC6,
        'SIC0':getSIC0,
        'SIC1':getSIC1,
        'SIC2':getSIC2,
        'SIC3':getSIC3,
        'SIC4':getSIC4,
        'SIC5':getSIC5,
        'SIC6':getSIC6,
        'IC0':getIC0,
        'IC1':getIC1,
        'IC2':getIC2,
        'IC3':getIC3,
        'IC4':getIC4,
        'IC5':getIC5,
        'IC6':getIC6
    }

def Getbasak(mol):
    """
    #################################################################
    Get the dictionary of basak descriptors for given moelcule mol
    #################################################################
    """
    dresult={}
    for DesLabel in _basak.keys():
        dresult[DesLabel]=round(_basak[DesLabel](mol),6)
    return dresult




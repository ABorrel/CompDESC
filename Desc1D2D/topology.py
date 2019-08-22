from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdchem
#from rdkit.Chem import pyPeriodicTable as PeriodicTable
from rdkit.Chem import GraphDescriptors as GD
import numpy
import scipy


periodicTable = rdchem.GetPeriodicTable()
print(periodicTable)

################################################################

def GetPrincipleQuantumNumber(atNum):
    """
    #################################################################
    *Internal Use Only*
    Get the principle quantum number of atom with atomic
    number equal to atNum 
    #################################################################
    """
    if atNum<=2:
        return 1
    elif atNum<=10:
        return 2
    elif atNum<=18:
        return 3
    elif atNum<=36:
        return 4
    elif atNum<=54:
        return 5
    elif atNum<=86:
        return 6
    else:
        return 7


def getWeiner(mol):
    return 1.0/2*sum(sum(Chem.GetDistanceMatrix(mol)))


def getMWeiner(mol):
    N = mol.GetNumAtoms()
    WeinerNumber = getWeiner(mol)
    return 2.0*WeinerNumber/(N*(N-1))


def getBalabanJ(mol):
    return Descriptors.BalabanJ(mol)


def getTigdi(mol):
    """
    #################################################################
    Calculation of graph distance index
    ---->Tigdi(log value)
    #################################################################
    """
    Distance= Chem.GetDistanceMatrix(mol)
    n = int(Distance.max())
    res = 0.0
    for i in range(n):
        temp=1./2*sum(sum(Distance==i+1))
        res = res+temp**2
    return res


def getLogTigdi(mol):
    """
    #################################################################
    Calculation of graph distance index
    ---->Tigdi(log value)
    #################################################################
    """
    Tigdi = getTigdi(mol)
    return numpy.log10(Tigdi)


def getXu(mol):
    nAT = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    Distance = Chem.GetDistanceMatrix(mol)
    sigma = scipy.sum(Distance,axis=1)
    temp1 = 0.0
    temp2 = 0.0
    for i in range(nAT):
        temp1=temp1+deltas[i]*((sigma[i])**2)
        temp2=temp2+deltas[i]*(sigma[i])
    Xu=numpy.sqrt(nAT)*numpy.log(temp1/temp2)
    return Xu





def CalculateDiameter(mol):
    Distance=Chem.GetDistanceMatrix(mol)
    return Distance.max()


def CalculateRadius(mol):
    Distance=Chem.GetDistanceMatrix(mol)
    temp=[]
    for i in Distance:
        temp.append(max(i))
    return min(temp)
    
def CalculatePetitjean(mol):
    """
    Value of (diameter - radius) / diameter as defined in [Petitjean 1992].
    ---->petitjeant
    """
    diameter=CalculateDiameter(mol)
    radius=CalculateRadius(mol)
    return 1-radius/float(diameter)


   


def getGMTI(mol):
    nAT = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    Distance = Chem.GetDistanceMatrix(mol)
    res=0.0
    for i in range(nAT):
        for j in range(i+1,nAT):
            res=res+deltas[i]*deltas[j]*Distance[i,j]
    return res

def getLogGMTI(mol):
    GMTI = getGMTI(mol)
    return numpy.log10(GMTI)
    
def getPol(mol):
    Distance= Chem.GetDistanceMatrix(mol)
    res = 1.0/2*sum(sum(Distance==3))
    return res


def getDZ(mol):
    """
    Calculation of Poglicani index
    The Pogliani index (Dz) is the sum over all non-hydrogen atoms
    of a modified vertex degree calculated as the ratio
    of the number of valence electrons over the principal
    quantum number of an atom [L. Pogliani, J.Phys.Chem.
    1996, 100, 18065-18077].
    """
    res = 0.0
    for atom in mol.GetAtoms():
        n = atom.GetAtomicNum()
        nV = periodicTable.GetNOuterElecs(n)
        mP = GetPrincipleQuantumNumber(n)
        res = res + ((nV)/mP)
    return res



def getIpc(mol):
    return Descriptors.Ipc(mol)

def getLogIpc(mol):
    ipc = getIpc(mol)
    return numpy.log10(ipc)

def CalculateBertzCT(mol):
    temp=GD.BertzCT(mol)
    if temp>0:
        return numpy.log10(temp)
    else:
        return "NaN"



def CalculateHarary(mol):

    Distance=numpy.array(Chem.GetDistanceMatrix(mol),'d')
                
    return 1.0/2*(sum(1.0/Distance[Distance!=0]))
        
    
def CalculateSchiultz(mol):
    Distance=numpy.array(Chem.GetDistanceMatrix(mol),'d')
    Adjacent=numpy.array(Chem.GetAdjacencyMatrix(mol),'d')
    VertexDegree=sum(Adjacent)
    
    return sum(scipy.dot((Distance+Adjacent),VertexDegree))



def CalculateZagreb1(mol):

    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    return sum(numpy.array(deltas)**2)


def CalculateZagreb2(mol):
    ke = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    
    return sum(ke)

def CalculateMZagreb1(mol):
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    res=sum((1./deltas)**2)
    return res
    

def CalculateMZagreb2(mol):
    cc = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc,'d')
    res = sum((1./cc)**2)
    return res

def CalculateQuadratic(mol):
    M=CalculateZagreb1(mol)
    N=mol.GetNumAtoms()
    return 3-2*N+M/2.0

def CalculatePlatt(mol):
    cc = [x.GetBeginAtom().GetDegree()+x.GetEndAtom().GetDegree()-2 for x in mol.GetBonds()]
    return sum(cc)



def CalculateSimpleTopoIndex(mol):
    """
    #################################################################
    Calculation of the logarithm of the simple topological index by Narumi,
    which is defined as the product of the vertex degree.
    ---->Sito
    #################################################################
    """
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    
    res=numpy.prod(deltas)
    if res>0:
        return numpy.log10(res)
    else:
        return "NaN"

def CalculateHarmonicTopoIndex(mol):
    """
    #################################################################
    Calculation of harmonic topological index proposed by Narnumi.
    ---->Hato
    #################################################################
    """
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')  
    nAtoms=mol.GetNumAtoms()
    
    res=nAtoms/sum(1./deltas)
    
    return res


def CalculateGeometricTopoIndex(mol):
    """
    #################################################################
    Geometric topological index by Narumi
    ---->Geto
    #################################################################
    """
    nAtoms=mol.GetNumAtoms()
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    
    temp=numpy.prod(deltas)
    res=numpy.power(temp,1./nAtoms)

    return res    

def CalculateArithmeticTopoIndex(mol):
    """
    #################################################################
    Arithmetic topological index by Narumi
    ---->Arto
    #################################################################
    """
    nAtoms=mol.GetNumAtoms()
    nBonds=mol.GetNumBonds()
    
    res=2.*nBonds/nAtoms
    return res



_topology={"Weiner": getWeiner,
           "Mweiner": getMWeiner,
           "BalabanJ": getBalabanJ,
           "Tigdi": getTigdi,
           "LogTigdi": getLogTigdi,
           "Xu": getXu,
           "GMTI": getGMTI,
           "LogGMTI": getLogGMTI,
           "Pol": getPol,
           "DZ": getDZ,
           "Ipc": getIpc}#,
           #"BertzCT": getBertzCT,
           #"Thara": getThara,
           #"Tsch": getTsch,
           #"ZM1": getZM1,
           #"ZM2": getZM2,
           #"MZM1": getMZM1,
           #"MZM2": getMZM2,
           #"Qindex": getQindex,
           #"Platt": getPlatt,
           #"diameter": getdiameter,
           #"radius": getradius,
           #"petitjean": getpetitjean,
           #"Sito": getSito,
           #"Hato": getHato,
           #"Geto": getGeto,
           #"Arto" :getArto}


def GetTopology(mol):
    dresult={}
    for DesLabel in _topology.keys():
        dresult[DesLabel] = round(_topology[DesLabel](mol), 6)
    return dresult



def _GetHTMLDoc():
    """
    #################################################################
    Write HTML documentation for this module.
    #################################################################
    """
    import pydoc
    pydoc.writedoc('topology')    
################################################################################

from rdkit import Chem
from rdkit.Chem import Descriptors, rdchem
import numpy
import scipy

periodicTable = rdchem.GetPeriodicTable()

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
    if Tigdi != 0.0:
        return numpy.log10(Tigdi)
    else:
        return 0.0


def getXu(mol):
    nAT = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    Distance = Chem.GetDistanceMatrix(mol)
    sigma = numpy.sum(Distance,axis=1)
    temp1 = 0.0
    temp2 = 0.0
    for i in range(nAT):
        temp1=temp1+deltas[i]*((sigma[i])**2)
        temp2=temp2+deltas[i]*(sigma[i])
    Xu=numpy.sqrt(nAT)*numpy.log(temp1/temp2)
    return Xu


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
    if GMTI != 0.0:
        return numpy.log10(GMTI)
    else:
        return 0.0
    
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
    if ipc != 0.0:
        return numpy.log10(ipc)
    else:
        return 0.0

def getBertzCT(mol):
    return Descriptors.BertzCT(mol)

def getLogBertzCT(mol):
    B = getBertzCT(mol)
    if B != 0.0:
        return numpy.log10(B)
    else:
        return 0.0

def getThara(mol):
    Distance = numpy.array(Chem.GetDistanceMatrix(mol),'d')
    return 1.0/2*(sum(1.0/Distance[Distance!=0]))
        
    
def getTsch(mol):
    Distance=numpy.array(Chem.GetDistanceMatrix(mol),'d')
    Adjacent=numpy.array(Chem.GetAdjacencyMatrix(mol),'d')
    VertexDegree=sum(Adjacent)
    return sum(numpy.dot((Distance+Adjacent),VertexDegree))

def getLogTsch(mol):
    T = getTsch(mol)
    if T != 0.0:
        return numpy.log10(T)
    else:
        return 0.0


def getZM1(mol):
    """
    based on the atom degree
    """
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    return sum(numpy.array(deltas)**2)


def getZM2(mol):
    """
    Based on the bond distance
    """
    ke = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    return sum(ke)

def getMZM1(mol):
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    res=sum((1./deltas)**2)
    return res
    

def getMZM2(mol):
    cc = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc,'d')
    res = sum((1./cc)**2)
    return res

def getQindex(mol):
    M=getZM1(mol)
    N=mol.GetNumAtoms()
    return 3-2*N+M/2.0

def getPlatt(mol):
    cc = [x.GetBeginAtom().GetDegree()+x.GetEndAtom().GetDegree()-2 for x in mol.GetBonds()]
    return sum(cc)


def getdiameterPJ(mol):
    Distance = Chem.GetDistanceMatrix(mol)
    return Distance.max()


def getradiusPJ(mol):
    Distance = Chem.GetDistanceMatrix(mol)
    temp = []
    for i in Distance:
        temp.append(max(i))
    return min(temp)


def getpetitjean(mol):
    """
    Value of (diameter - radius) / diameter as defined in [Petitjean 1992].
    ---->petitjeant
    """
    diameter = getdiameterPJ(mol)
    radius = getradiusPJ(mol)
    return 1.0 - radius / float(diameter)


def getSito(mol):
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
    
    res = numpy.prod(deltas)
    if res > 0:
        return numpy.log10(res)
    else:
        return "NA"

def getHato(mol):
    """
    #################################################################
    Calculation of harmonic topological index proposed by Narnumi.
    ---->Hato
    #################################################################
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0.0 in deltas:
        deltas.remove(0.0)
    deltas = numpy.array(deltas,'d')
    nAtoms = mol.GetNumAtoms()

    den = sum(1./deltas)
    if den == 0.0:
        res =  0.0
    else:
        res = nAtoms/den
    return res


def getGeto(mol):
    """
    #################################################################
    Geometric topological index by Narumi
    ---->Geto
    #################################################################
    """
    nAtoms = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas,'d')
    
    temp = numpy.prod(deltas)
    res = numpy.power(temp,1./nAtoms)

    return res    

def getArto(mol):
    """
    #################################################################
    Arithmetic topological index by Narumi
    ---->Arto
    #################################################################
    """
    nAtoms = mol.GetNumAtoms()
    nBonds = mol.GetNumBonds()
    
    res = 2.0*nBonds/nAtoms
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
           "Ipc": getIpc,
           "LogIpc": getLogIpc,
           "BertzCT": getBertzCT,
           "LogBertzCT": getLogBertzCT,
           "Thara": getThara,
           "Tsch": getTsch,
           "LogTsch": getLogTsch,
           "ZM1": getZM1,
           "ZM2": getZM2,
           "MZM1": getMZM1,
           "MZM2": getMZM2,
           "Qindex": getQindex,
           "Platt": getPlatt,
           "diameterPJ": getdiameterPJ,
           "radiusPJ": getradiusPJ,
           "petitjean": getpetitjean,
           "Sito": getSito,
           "Hato": getHato,
           "Geto": getGeto,
           "Arto" :getArto}


def GetTopology(mol):
    dresult={}
    for DesLabel in _topology.keys():
        dresult[DesLabel] = round(_topology[DesLabel](mol), 6)
    return dresult



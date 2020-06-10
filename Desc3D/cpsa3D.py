from .asa3D import calculate_asa
from .Atom3DProperty import GetAtomClassList

def GetChargeSA(lcoordinates, RadiusProbe=1.5, n_sphere_point=960):
    """Get the list form for all atoms in a molecule.
    It includes the atom symbol, charge and partial solvent-accessible
    surface areas.
    Note that this is list form whose element is still list form of each atom.

    """
    atoms=GetAtomClassList(lcoordinates)
    FASA=calculate_asa(atoms, RadiusProbe, n_sphere_point)
    res=[]
    for i in range(len(FASA)):
        res.append([lcoordinates[i][3],lcoordinates[i][4],FASA[i]])
    return res

#####################################

def getASA(ChargeSA):
    """ The calculation of solvent-accessible surface areas
    -->ASA """
    res=0.0
    for i in ChargeSA:
        res = res + i[2]
    return res


def getMSA(lcoordinates):
    """ The calculation of molecular surface areas
    -->MSA"""
    ChargeSA=GetChargeSA(lcoordinates, RadiusProbe=0, n_sphere_point=960)
    res=0.0
    for i in ChargeSA:
        res=res+i[2]
    return res


def getPNSA1(ChargeSA):
    """ The calculation of partial negative area
    It is the sum of the solvent-accessible surface areas of all 
    negatively charged atoms.
    -->PNSA1"""
    res=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            res=res+i[2]

    return res


def getPPSA1(ChargeSA):
    """ The calculation of partial negative area
    It is the sum of the solvent-accessible surface areas of
    all positively charged atoms.
    -->PPSA1
    """
    res=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            res=res+i[2]

    return res

def getPNSA2(ChargeSA):
    """The calculation of total charge wighted negative surface area
    It is the partial negative solvent-accessible surface area
    multiplied by the total negative charge.
    -->PNSA2
    """
    temp1=temp2=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            temp1=temp1+float(i[1])
            temp2=temp2+i[2]
    res=temp1*temp2

    return res


def getPPSA2(ChargeSA):
    """The calculation of total charge wighted negative surface area
    It is the partial negative solvent-accessible surface area 
    multiplied by the total positive charge.
    -->PPSA2 """

    temp1=temp2=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            temp1=temp1+float(i[1])
            temp2=temp2+i[2]
    res=temp1*temp2

    return res

def getPNSA3(ChargeSA):
    """The calculation of atom charge weighted negative surface ares
    It is the sum of the products of atomic solvent-accessible 
    surface area and partial charges over all negatively charges atoms.
    -->PNSA3 """

    res=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            res=res+float(i[1])*i[2]
    return res


def getPPSA3(ChargeSA):
    """The calculation of atom charge weighted positive surface ares
    It is the sum of the products of atomic solvent-accessible
    surface area and partial charges over all positively charges atoms.
    -->PPSA3"""

    res=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            res=res+float(i[1])*i[2]
    return res


def getDPSA1(ChargeSA):
    """ The calculation of difference in charged partial surface area
    -->DPSA1"""
    return getPPSA1(ChargeSA)-getPNSA1(ChargeSA)

def getDPSA2(ChargeSA):
    """ The calculation of difference in charged partial surface area
    -->DPSA2"""
    return getPPSA2(ChargeSA)-getPNSA2(ChargeSA)

def getDPSA3(ChargeSA):
    """ The calculation of difference in charged partial surface area
    -->DPSA3"""
    return getPPSA3(ChargeSA)-getPNSA3(ChargeSA)

def getFNSA1(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FNSA1"""
    temp = 0.0
    for i in ChargeSA:
        temp = temp + i[2]
    if temp == 0.0:
        return 0.0
    return getPNSA1(ChargeSA)/temp

def getFNSA2(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FNSA2"""
    temp = 0.0
    for i in ChargeSA:
        temp = temp + i[2]
    if temp == 0.0:
        return 0.0
    return getPNSA2(ChargeSA)/temp

def getFNSA3(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FNSA3"""
    temp = 0.0
    for i in ChargeSA:
        temp = temp + i[2]
    if temp == 0.0:
        return 0.0
    return getPNSA3(ChargeSA)/temp


def getFPSA1(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FPSA1
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    if temp == 0.0:
        return 0.0
    else:
        return getPPSA1(ChargeSA)/temp

def getFPSA2(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FPSA2
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    if temp == 0.0:
        return 0.0
    else:
        return getPPSA2(ChargeSA)/temp

def getFPSA3(ChargeSA):
    """The calculation of fractional charged partial negative surface areas
    -->FPSA3
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    if temp == 0.0:
        return 0.0
    else:
        return getPPSA3(ChargeSA)/temp


def getWNSA1(ChargeSA):
    """The calculation of surface weighted charged partial negative 
    surface areas
    -->WNSA1
    """
    temp = 0.0
    for i in ChargeSA:
        temp = temp+i[2]
    if temp == 0.0:
        return 0.0

    return getPNSA1(ChargeSA)*temp/1000

def getWNSA2(ChargeSA):
    """The calculation of surface weighted charged partial negative
    surface areas
    -->WNSA2
    """
    temp = 0.0
    for i in ChargeSA:
        temp = temp+i[2]
    if temp == 0.0:
        return 0.0

    return getPNSA2(ChargeSA)*temp/1000.0

def getWNSA3(ChargeSA):
    """The calculation of surface weighted charged partial negative
    surface areas
    -->WNSA3
    """
    temp = 0.0
    for i in ChargeSA:
        temp = temp+i[2]
    if temp == 0.0:
        return 0.0

    return getPNSA3(ChargeSA)*temp/1000.0


def getWPSA1(ChargeSA):
    """The calculation of surface weighted charged partial negative 
    surface areas
    -->WPSA1
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    if temp == 0.0:
        return 0.0
    return getPPSA1(ChargeSA)*temp/1000.0

def getWPSA2(ChargeSA):
    """The calculation of surface weighted charged partial negative
    surface areas
    -->WPSA2
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    if temp == 0.0:
        return 0.0
    return getPPSA2(ChargeSA)*temp/1000.0

def getWPSA3(ChargeSA):
    """The calculation of surface weighted charged partial negative
    surface areas
    -->WPSA3
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    if temp == 0.0:
        return 0.0
    return getPPSA3(ChargeSA)*temp/1000.0

def getTASA(ChargeSA):
    """The calculation of total hydrophobic surface area
    -->TASA
    """
    res = 0.0
    for i in ChargeSA:
        if abs(float(i[1])) < 0.2:
            res = res + i[2]
    return res


def getPSA(ChargeSA):
    """The calculation of total polar surface area
    -->PSA
    """
    res=0.0
    for i in ChargeSA:
        if abs(float(i[1]))>=0.2:
            res=res+i[2]
    return res


def getTATP(ChargeSA):
    """The fraction between TASA and TPSA
    --->FrTATP
    """
    res=0.0
    if getPSA(ChargeSA) == 0.0:
        return res
    else:
        return getTASA(ChargeSA)/getPSA(ChargeSA)


def getRASA(ChargeSA):
    """The calculation of relative hydrophobic surface area
    -->RASA
    """
    temp=0.0
    for i in ChargeSA:
        temp = temp + i[2]
    if temp == 0.0:
        return 0.0
    return getTASA(ChargeSA)/temp


def getRPSA(ChargeSA):
    """The calculation of relative polar surface area
    -->RPSA
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    if temp == 0.0:
        return 0.0
    return getPSA(ChargeSA)/temp


def getRNCS(ChargeSA):
    """The calculation of relative negative charge surface area
    -->RNCS
    """
    charge=[]
    for i in ChargeSA:
        charge.append(float(i[1]))

    temp=[]
    for i in ChargeSA:
        temp.append(i[2])

    try:
        RNCG = min(charge)/sum([i for i in charge if i < 0.0])
        return temp[charge.index(min(charge))]/RNCG
    except:
        return 0.0


def getRPCS(ChargeSA):
    """The calculation of relative positive charge surface area
    -->RPCS
    """
    charge=[]
    for i in ChargeSA:
        charge.append(float(i[1]))

    temp=[]
    for i in ChargeSA:
        temp.append(i[2])

    try:
        RPCG=max(charge)/sum([i for i in charge if i>0])
        return temp[charge.index(min(charge))]/RPCG
    except:
        return 0.0


_cpsa3D = {"ASA": getASA,
           "MSA": getMSA,
           "PNSA1": getPNSA1,
           "PNSA2": getPNSA2,
           "PNSA3": getPNSA3,
           "PPSA1": getPPSA1,
           "PPSA2": getPPSA2,
           "PPSA3": getPPSA3,
           "DPSA1": getDPSA1,
           "DPSA2": getDPSA2,
           "DPSA3": getDPSA3,
           "FNSA1": getFNSA1,
           "FNSA2": getFNSA2,
           "FNSA3": getFNSA3,
           "FPSA1": getFPSA1,
           "FPSA2": getFPSA2,
           "FPSA3": getFPSA3,
           "WNSA1": getWNSA1,
           "WNSA2": getWNSA2,
           "WNSA3": getWNSA3,
           "WPSA1": getWPSA1,
           "WPSA2": getWPSA2,
           "WPSA3": getWPSA3,
           "TASA": getTASA,
           "PSA": getPSA,
           "RASA": getRASA,
           "RPSA": getRPSA,
           "RNCS": getRNCS,
           "RPCS": getRPCS,
           "TATP": getTATP}


import math
def GetCPSA3D(lcoordinates):
    ChargeSA = GetChargeSA(lcoordinates, RadiusProbe=1.5, n_sphere_point=5000)
    dresult = {}
    for DesLabel in _cpsa3D.keys():
        if DesLabel == "MSA":
            val = _cpsa3D[DesLabel](lcoordinates)
        else:
            val =_cpsa3D[DesLabel](ChargeSA)

        if math.isnan(val) == True:
            dresult[DesLabel] = "NA"
        else:
            dresult[DesLabel] = round(val, 6)

    return dresult


#########################
# GEOMETRY DESCRIPTORS  #
#########################
import scipy
import scipy.linalg
import math
from .Atom3DProperty import get_atomicMass

from rdkit.Chem import Descriptors3D
import numpy as np

############################################################################


def GetAtomDistance(x, y):
    """
    #################################################################
    Obtain the Elucidian distance based on the coordinates of two atoms
    #################################################################
    """

    temp = [math.pow(x[0] - y[0], 2), math.pow(x[1] - y[1], 2), math.pow(x[2] - y[2], 2)]
    res = np.sqrt(sum(temp))
    return res


def GetGementricalDistanceMatrix(CoordinateList):
    """
    #################################################################
    Obtain the distance matrix of a molecule based on coordinate list
    #################################################################
    """
    NAtom = len(CoordinateList)
    DistanceMatrix = np.zeros((NAtom, NAtom))
    for i in range(NAtom - 1):
        for j in range(i + 1, NAtom):
            DistanceMatrix[i, j] = GetAtomDistance(CoordinateList[i], CoordinateList[j])
            DistanceMatrix[j, i] = DistanceMatrix[i, j]
    return DistanceMatrix


def _GetMassCenter(MassCoordinates):
    """
    #################################################################
    Get the center of mass.
    INPUT: MassCoordinates is [[atommass,[x,y,z]],......].
    #################################################################
    """

    res1 = 0.0
    res2 = 0.0
    res3 = 0.0
    temp = []
    for i in MassCoordinates:
        res1 = res1 + i[0] * i[1][0]
        res2 = res2 + i[0] * i[1][1]
        res3 = res3 + i[0] * i[1][2]
        temp.append(i[0])
    result = [res1 / sum(temp), res2 / sum(temp), res3 / sum(temp)]
    return result


def _GetGeometricalCenter(ChargeCoordinates):
    """
    #################################################################
    Get the geometrical center
    #################################################################
    """
    res1 = []
    res2 = []
    res3 = []

    for i in ChargeCoordinates:
        res1.append(float(i[0]))
        res2.append(float(i[1]))
        res3.append(float(i[2]))

    result = [scipy.mean(res1), scipy.mean(res2), scipy.mean(res3)]

    return result

def GetInertiaMatrix(ChargeCoordinates):
    """
    #################################################################
    Get Inertia matrix based on atomic mass and optimized coordinates.
    #################################################################
    """
    temp = []
    for coords in ChargeCoordinates:
        temp.append([get_atomicMass(coords[3]), [coords[0], coords[1], coords[2]]])


    nAT = len(temp)
    InertiaMatrix = scipy.zeros((3, 3))
    res11 = 0.0
    res22 = 0.0
    res33 = 0.0
    res12 = 0.0
    res23 = 0.0
    res13 = 0.0
    for i in range(nAT):
        res11 = res11 + temp[i][0] * (math.pow(temp[i][1][1], 2) + math.pow(temp[i][1][2], 2))
        res22 = res22 + temp[i][0] * (math.pow(temp[i][1][0], 2) + math.pow(temp[i][1][2], 2))
        res33 = res33 + temp[i][0] * (math.pow(temp[i][1][0], 2) + math.pow(temp[i][1][1], 2))
        res12 = res12 + temp[i][0] * (temp[i][1][0] * temp[i][1][1])
        res13 = res13 + temp[i][0] * (temp[i][1][0] * temp[i][1][2])
        res23 = res23 + temp[i][0] * (temp[i][1][1] * temp[i][1][2])
    InertiaMatrix[0, 0] = res11
    InertiaMatrix[1, 1] = res22
    InertiaMatrix[2, 2] = res33
    InertiaMatrix[0, 1] = res12
    InertiaMatrix[0, 2] = res13
    InertiaMatrix[1, 2] = res23
    InertiaMatrix[1, 0] = res12
    InertiaMatrix[2, 0] = res13
    InertiaMatrix[2, 1] = res23

    return InertiaMatrix


def CalculatePrincipalMomentofInertia(mol, ChargeCoordinates):
    """
    #################################################################
    X,Y and Z-principal geometric moment.
    drived from ADAPT developed by Jurs.
    #################################################################
    """
    InertiaMatrix = GetInertiaMatrix(ChargeCoordinates)
    ma = scipy.mean(InertiaMatrix, axis=1)
    ms = scipy.std(InertiaMatrix, axis=1, ddof=1)
    bb = scipy.ones((3, 1))
    InertiaMatrix = (InertiaMatrix - bb * ma.T) / (bb * ms.T)
    u, s, v = scipy.linalg.svd(InertiaMatrix)

    res = {}
    res['IA'] = round(s[2], 3)
    res['IB'] = round(s[1], 3)
    res['IC'] = round(s[0], 3)

    return res


def CalculateRatioPMI(mol, ChargeCoordinates):
    """
    #################################################################
    The ratio of X/Y, Y/Z and X/Z (principal moment of inertia)
    drived from ADAPT developed by Jurs.
    #################################################################
    """
    temp = CalculatePrincipalMomentofInertia(mol, ChargeCoordinates)
    res = {}

    res['IA/B'] = round(temp['IA'] / temp['IB'], 3)
    res['IA/C'] = round(temp['IA'] / temp['IC'], 3)
    res['IB/C'] = round(temp['IB'] / temp['IC'], 3)
    return res


################
### for desc ###
################

def getW3DH(coords):
    """
    #################################################################
    The calculation of 3-D Wiener index based geometrical distance matrix optimized
    by MOPAC(including Hs)
    -->W3DH
    #################################################################
    """
    temp = []
    for i in coords:
        temp.append([float(i[0]), float(i[1]), float(i[2])])
    DistanceMatrix = GetGementricalDistanceMatrix(temp)
    return np.sum(DistanceMatrix) / 2.0


def getW3D(coords):
    """
    #################################################################
    The calculation of 3-D Wiener index based
    gemetrical distance matrix optimized
    by MOPAC(Not including Hs)
    -->W3D
    #################################################################
    """
    temp = []
    for i in coords:
        if i[3] != 'H':
            temp.append([float(i[0]), float(i[1]), float(i[2])])
    DistanceMatrix = GetGementricalDistanceMatrix(temp)
    return np.sum(DistanceMatrix) / 2.0


def getPetitj3D(ChargeCoordinates):
    """
    #################################################################
    CalculatePetitjean Index based on molecular gemetrical distance matrix
    -->Petitj3D
    The 3D Petitjean shape index (PJI3) is calculated
    dividing the difference between geometric diameter and
    radius by the geometric radius [P.A. Bath, A.R. Poirrette,
    P. Willett, F.H. Allen, J.Chem.Inf.Comput.Sci. 1995, 35, 714-716].
    The geometric radius of a molecule is defined as the minimum
    geometric eccentricity and the diameter is defined as the
    maximum geometric eccentricity in the molecule, the atom
    geometric eccentricity being the longest geometric distance
    from the considered atom to any other atom in the molecule.
    #################################################################
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[0]), float(i[1]), float(i[2])])

    DistanceMatrix = GetGementricalDistanceMatrix(temp)
    temp1 = np.amax(DistanceMatrix, axis=0)

    return max(temp1) / min(temp1) - 1.0


def getGeDi(ChargeCoordinates):
    """
    #################################################################
    The longest distance between two atoms (gemetrical diameter)
    -->GeDi
    #################################################################
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[0]), float(i[1]), float(i[2])])
    DistanceMatrix = GetGementricalDistanceMatrix(temp)
    temp1 = np.amax(DistanceMatrix, axis=0)
    return max(temp1)


def getgrav(ChargeCoordinates):
    """
    #################################################################
    Calculation of Gravitational 3D index.
    --->grav
    #################################################################
    """
    temp = []
    for coords in ChargeCoordinates:
        temp.append([get_atomicMass(coords[3]), [coords[0], coords[1], coords[2]]])

    nAT = len(temp)
    result = 0.0
    for i in range(nAT - 1):
        for j in range(i + 1, nAT):
            dis = GetAtomDistance(temp[i][1], temp[j][1])
            den = np.power(dis, 2)
            if den != 0.0:
                result = result + temp[i][0] * temp[j][0] / np.power(dis, 2)
    return float(result) / 100.0


def getRadiusOfGyration(mol3D):
    """
    #################################################################
    Calculation of Radius of gyration.
    --->rygr
    #################################################################
    """
    # case of no conformer generated
    try: return Descriptors3D.RadiusOfGyration(mol3D)
    except: return 0.0

def getHarary3D(ChargeCoordinates):
    """
    #################################################################
    The 3D-Harary index (H3D) is calculated as
    the sum of all the reciprocal geometric distances
    in a molecule.
    --->Harary3D
    #################################################################
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[0]), float(i[1]), float(i[2])])
    DistanceMatrix = GetGementricalDistanceMatrix(temp)
    nAT = len(temp)
    res = 0.0
    for i in range(nAT - 1):
        for j in range(i + 1, nAT):
            if DistanceMatrix[i, j] == 0:
                cds = 0.0
            else:
                cds = 1. / DistanceMatrix[i, j]
            res = res + cds
    return res


def getAGDD(ChargeCoordinates):
    """
    #################################################################
    The average geometric distance degree (AGDD) is
    calculated dividing the sum of all geometric distance
    degrees by the total number of molecule atoms (nAT).
    ---->AGDD
    #################################################################
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[0]), float(i[1]), float(i[2])])
    DistanceMatrix = GetGementricalDistanceMatrix(temp)
    nAT = len(temp)
    res = sum(sum(DistanceMatrix)) / nAT
    return res


def getSEig(ChargeCoordinates):
    """
    #################################################################
    The absolute eigenvalue sum on geometry matrix (SEig)
    is the sum of the absolute eigenvalues of the geometry matrix.
    --->SEig
    #################################################################
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[0]), float(i[1]), float(i[2])])
    DistanceMatrix = GetGementricalDistanceMatrix(temp)
    u, s, vt = scipy.linalg.svd(DistanceMatrix)
    return sum(abs(s))


def getSPAN(ChargeCoordinates):
    """
    #################################################################
    The span R (SPAN) is a size descriptor defined as
    the radius of the smallest sphere, centred on the centre
    of mass, completely enclosing all atoms of a molecule
    [G.A. Arteca, Molecular Shape Descriptors in Reviews in
    Computational Chemistry - Vol. 9, K.B. Lipkowitz, D. Boyd (Eds.),
    VCH Publishers, New York (NY), pp. 191-253, 1991]
    --->SPAN
    #################################################################
    """
    temp = []
    for coords in ChargeCoordinates:
        temp.append([get_atomicMass(coords[3]), [coords[0], coords[1], coords[2]]])
    masscenter = _GetMassCenter(temp)
    res = []
    for i in temp:
        res.append(GetAtomDistance(i[1], masscenter))
    return float(max(res))


def getASPAN(ChargeCoordinates):
    """
    #################################################################
    The average span R (SPAM) is the root square of
    the ratio of SPAN over the number of atoms.
    --->ASPAN
    #################################################################
    """
    temp = []
    for coords in ChargeCoordinates:
        temp.append([get_atomicMass(coords[3]), [coords[0], coords[1], coords[2]]])
    nAT = len(temp)
    masscenter = _GetMassCenter(temp)
    res = []
    for i in temp:
        res.append(GetAtomDistance(i[1], masscenter))
    return math.pow(float(max(res)) / nAT, 0.5)


def getEccentricity(mol3D):
    """
    #################################################################
    The molecular eccentricity (MEcc) is a shape descriptor
    calculated from the eigenvalues l of the molecular inertia matrix
    [G.A. Arteca, Molecular Shape Descriptors in Reviews
    in Computational Chemistry - Vol. 9, K.B. Lipkowitz, D. Boyd (Eds.),
    VCH Publishers, New York (NY), pp. 191-253, 1991].
    --->MEcc
    #################################################################
    """
    try: return Descriptors3D.Eccentricity(mol3D)
    except: return "NA"

def getAsphericity(mol3D):
    try: return Descriptors3D.Asphericity(mol3D)
    except: return "NA"

def getInertialShapeFactor(mol3D):
    try: return Descriptors3D.InertialShapeFactor(mol3D)
    except: return "NA"

def getNPR1(mol3D):
    try: return Descriptors3D.NPR1(mol3D)
    except: return "NA"

def getNPR2(mol3D):
    try: return Descriptors3D.NPR2(mol3D)
    except: return "NA"

def getPMI1(mol3D):
    try: return Descriptors3D.PMI1(mol3D)
    except: return "NA"

def getPMI2(mol3D):
    try: return Descriptors3D.PMI2(mol3D)
    except: return "NA"

def getPMI3(mol3D):
    try: return Descriptors3D.PMI3(mol3D)
    except: return "NA"

def getSpherocityIndex(mol3D):
    try: return Descriptors3D.SpherocityIndex(mol3D)
    except: return "NA"


_geo3D = {"W3DH": getW3DH,
          "W3D": getW3D,
          "Petitj3D": getPetitj3D,
          "GeDi": getGeDi,
          "grav": getgrav,
          "RadiusOfGyration": getRadiusOfGyration,
          "Harary3D": getHarary3D,
          "AGDD": getAGDD,
          "SEig": getSEig,
          "SPAN": getSPAN,
          "ASPAN": getASPAN,
          "Eccentricity": getEccentricity,
          "Asphericity": getAsphericity,
          "InertialShapeFactor": getInertialShapeFactor,
          "NPR1": getNPR1,
          "NPR2": getNPR2,
          "PMI1": getPMI1,
          "PMI2": getPMI2,
          "PMI3": getPMI3,
          "SpherocityIndex": getSpherocityIndex}



def GetGeo3D(coords, mol3D):
    dresult = {}
    for DesLabel in _geo3D.keys():
        if DesLabel == "W3DH" or DesLabel == "W3D" or DesLabel == "Petitj3D" or DesLabel == "GeDi" or DesLabel == "grav" or DesLabel == "Harary3D" or DesLabel == "AGDD" or DesLabel == "SEig" or DesLabel == "SPAN" or DesLabel == "ASPAN":
            val = _geo3D[DesLabel](coords)
        else:
            val = _geo3D[DesLabel](mol3D)

        if val == "NA":
            dresult[DesLabel] = "NA"
        elif math.isnan(val) == True:
            dresult[DesLabel] = "NA"
        else:
            dresult[DesLabel] = round(val, 6)
    return dresult


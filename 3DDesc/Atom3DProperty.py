import scipy
from re import search
from .vector3d import Vector3d

class Atom:
    """
    #################################################################
    A atom class used for wrapping some properties of atoms.

    Note that Coordinates is the output of the function

    (_ReadCoordinates).
    #################################################################
    """

    def __init__(self, Coordinates):

        self.pos = Vector3d()
        self.radius = 0.0
        self.Coordinates = Coordinates
        self.Element = ''

    def SetCoordinates(self):

        temp = self.Coordinates
        self.pos.x = float(temp[0])
        self.pos.y = float(temp[1])
        self.pos.z = float(temp[2])

    def GetCoordinates(self):

        self.SetCoordinates()

        return self.pos

    def SetElement(self):

        temp = self.Coordinates

        self.Element = temp[3]

    def GetElement(self):

        self.SetElement()

        return self.Element

    def SetRadius(self):

        radii = {'H': 1.20, 'N': 1.55, 'Na': 2.27, 'Cu': 1.40, 'Cl': 1.75, 'C': 1.70,
                 'O': 1.52, 'I': 1.98, 'P': 1.80, 'B': 1.85, 'Br': 1.85, 'S': 1.80, 'Se': 1.90,
                 'F': 1.47, 'Fe': 1.80, 'K': 2.75, 'Mn': 1.73, 'Mg': 1.73, 'Zn': 1.39, 'Hg': 1.8,
                 'Li': 1.8, '.': 1.8}

        temp = self.GetElement()

        if temp in radii.keys():
            self.radius = radii[temp]
        else:
            self.radius = radii['.']

    def GetRadius(self):

        self.SetRadius()

        return self.radius




###########################################################################

def GetAtomClassList(Coordinates):
    """
    #################################################################
    Combine all atoms in a molecule into a list form.
    Note that Coordinates is the output of the function (_ReadCoordinates).
    #################################################################
    """
    Atoms = []
    for i in Coordinates:
        atom = Atom(i)
        atom.SetCoordinates()
        atom.SetElement()
        atom.SetRadius()
        Atoms.append(atom)
    return Atoms


def GetAtomCoordinateMatrix(lcoordinates):
    """
    #################################################################
    Get the atom coordinate matrix
    #################################################################
    """
    nAtom = len(lcoordinates)
    CoordinateMatrix = scipy.zeros([nAtom, 3])
    AtomLabel = []

    for i, j in enumerate(lcoordinates):
        CoordinateMatrix[i, :] = [j[0], j[1], j[2]]
        AtomLabel.append(j[3])

    return scipy.matrix(CoordinateMatrix), AtomLabel



def get_atomicMass(element):

    atomicMass={'H': 1.0079, 'N': 14.0067, 'Na': 22.9897, 'Cu': 63.546, 'Cl': 35.453, 'C': 12.0107,
                 'O': 15.9994, 'I': 126.9045, 'P': 30.9738, 'B': 10.811, 'Br': 79.904, 'S': 32.065, 'Se': 78.96,
                 'F': 18.9984, 'Fe': 55.845, 'K': 39.0983, 'Mn': 54.938, 'Mg': 24.305, 'Zn': 65.39, 'Hg': 200.59,
                 'Li': 6.941, 'Co': 58.9332, "Si":28.0855, "As": 74.9216, "Te":127.6, "Sr":87.62}


    return atomicMass[element]


def get_MW(lcoords, H=1):
    MW = 0
    for coords in lcoords:
        if coords[3] != "H" and H ==0:MW = MW + get_atomicMass(coords[3])
        else:MW = MW + get_atomicMass(coords[3])
    return MW

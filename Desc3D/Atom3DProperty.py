import scipy
from .vector3d import Vector3d
from rdkit.Chem import rdchem
periodicTable = rdchem.GetPeriodicTable()

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

        temp = self.GetElement()
        self.radius = periodicTable.GetRvdw(temp)
        #except:self.radius = 1.8

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
    return periodicTable.GetAtomicWeight(element)


def get_MW(lcoords, H=1):
    MW = 0
    for coords in lcoords:
        if coords[3] != "H" and H ==0:MW = MW + get_atomicMass(coords[3])
        else:MW = MW + get_atomicMass(coords[3])
    return MW

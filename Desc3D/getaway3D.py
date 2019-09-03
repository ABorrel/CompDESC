from rdkit.Chem import rdMolDescriptors

def getGETAWAY(mol3D):
    dout = {}
    lgetaway = rdMolDescriptors.CalcGETAWAY(mol3D)
    for i in range(1, len(lgetaway) + 1):
        dout["GETAWAY" + str(i)] = round(lgetaway[i - 1], 6)
    return dout


_getaway3D = {}
ldesc = ["GETAWAY" + str(i) for i in range(1, 274)]
for desc in ldesc:
    _getaway3D[desc] = getGETAWAY


def GetGETAWAY3D(mol3D):
    return getGETAWAY(mol3D)
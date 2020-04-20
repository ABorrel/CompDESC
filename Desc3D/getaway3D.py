from rdkit.Chem import rdMolDescriptors
import math
import numpy

def getGETAWAY(mol3D):
    dout = {}

    try: lgetaway = rdMolDescriptors.CalcGETAWAY(mol3D)
    except: lgetaway = []
    if lgetaway != []:
        for i in range(1, len(lgetaway) + 1):
            val = lgetaway[i - 1]
            if math.isnan(val):
                dout["GETAWAY" + str(i)] = 0.0
            elif val == numpy.inf:
                dout["GETAWAY" + str(i)] = "NA"
            else:
                dout["GETAWAY" + str(i)] = round(val, 6)
        return dout
    else:
        for desc in _getaway3D.keys():
            dout[desc] = "NA"
        return dout



_getaway3D = {}
ldesc = ["GETAWAY" + str(i) for i in range(1, 274)]
for desc in ldesc:
    _getaway3D[desc] = getGETAWAY


def GetGETAWAY3D(mol3D):
    return getGETAWAY(mol3D)
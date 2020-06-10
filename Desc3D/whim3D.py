from rdkit.Chem import rdMolDescriptors
import math

def getWHIM(mol3D):
    dout = {}
    try: lwhim = rdMolDescriptors.CalcWHIM(mol3D)
    except: lwhim = []

    if lwhim != []:
        for i in range(1, len(lwhim) + 1):
            val = lwhim[i - 1]
            if math.isnan(val) == True:
                dout["WHIM" + str(i)] = "NA"
            else:
                dout["WHIM" + str(i)] = round(val, 6)
        return dout
    else:
        for desc in _whim3D.keys():
            dout[desc] = "NA"
        return dout


_whim3D = {}
ldesc = ["WHIM" + str(i) for i in range(1, 115)]
for desc in ldesc:
    _whim3D[desc] = getWHIM


def GetWHIM3D(mol3D):
    return getWHIM(mol3D)
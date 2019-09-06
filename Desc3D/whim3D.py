from rdkit.Chem import rdMolDescriptors

def getWHIM(mol3D):
    dout = {}
    lwhim = rdMolDescriptors.CalcWHIM(mol3D)

    try: lwhim = rdMolDescriptors.CalcWHIM(mol3D)
    except: lwhim = []
    if lwhim != []:
        for i in range(1, len(lwhim) + 1):
            dout["WHIM" + str(i)] = round(lwhim[i - 1], 6)
        return dout
    else:
        for desc in _whim3D:
            dout[desc] = "NA"
        return dout


_whim3D = {}
ldesc = ["WHIM" + str(i) for i in range(1, 115)]
for desc in ldesc:
    _whim3D[desc] = getWHIM


def GetWHIM3D(mol3D):
    return getWHIM(mol3D)
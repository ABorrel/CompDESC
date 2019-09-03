from rdkit.Chem import rdMolDescriptors

def getAUTOCORR(mol3D):
    dout = {}
    lautocorr3D = rdMolDescriptors.CalcAUTOCORR3D(mol3D)
    for i in range(1, len(lautocorr3D) + 1):
        dout["AUTOCORR3D" + str(i)] = round(lautocorr3D[i - 1], 6)
    return dout


_autocorr3D = {}
ldesc = ["AUTOCORR3D" + str(i) for i in range(1, 81)]
for desc in ldesc:
    _autocorr3D[desc] = getAUTOCORR


def GetAUTOCORR3D(mol3D):
    return getAUTOCORR(mol3D)
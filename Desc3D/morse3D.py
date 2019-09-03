from rdkit.Chem import rdMolDescriptors

def getMORSE(mol3D):
    dout = {}
    lmorse = rdMolDescriptors.CalcMORSE(mol3D)
    for i in range(1, len(lmorse) + 1):
        dout["MORSE" + str(i)] = round(lmorse[i - 1], 6)
    return dout


_morse3D = {}
ldesc = ["MORSE" + str(i) for i in range(1, 225)]
for desc in ldesc:
    _morse3D[desc] = getMORSE


def GetMorse3D(mol3D):
    return getMORSE(mol3D)

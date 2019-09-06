from rdkit.Chem import rdMolDescriptors

def getMORSE(mol3D):
    dout = {}
    try: lmorse = rdMolDescriptors.CalcMORSE(mol3D)
    except: lmorse = []
    if lmorse != []:
        for i in range(1, len(lmorse) + 1):
            dout["MORSE" + str(i)] = round(lmorse[i - 1], 6)
        return dout
    else:
        for desc in _morse3D:
            dout[desc] = "NA"
        return dout


_morse3D = {}
ldesc = ["MORSE" + str(i) for i in range(1, 225)]
for desc in ldesc:
    _morse3D[desc] = getMORSE


def GetMorse3D(mol3D):
    return getMORSE(mol3D)

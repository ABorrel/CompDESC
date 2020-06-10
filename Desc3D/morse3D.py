from rdkit.Chem import rdMolDescriptors
import math

def getMORSE(mol3D):
    dout = {}
    try: lmorse = rdMolDescriptors.CalcMORSE(mol3D)
    except: lmorse = []

    if lmorse != []:
        for i in range(1, len(lmorse) + 1):
            val = lmorse[i - 1]
            if math.isnan(val) == True:
                dout["MORSE" + str(i)] = "NA"
            else:
                dout["MORSE" + str(i)] = round(val, 6)

        return dout
    else:
        for desc in _morse3D.keys():
            dout[desc] = "NA"
        return dout


_morse3D = {}
ldesc = ["MORSE" + str(i) for i in range(1, 225)]
for desc in ldesc:
    _morse3D[desc] = getMORSE


def GetMorse3D(mol3D):
    return getMORSE(mol3D)

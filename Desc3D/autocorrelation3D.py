from rdkit.Chem import rdMolDescriptors
import math

def getAUTOCORR(mol3D):
    dout = {}
    try: lautocorr3D = rdMolDescriptors.CalcMORSE(mol3D)
    except: lautocorr3D = []

    if lautocorr3D != []:
        for i in range(1, len(lautocorr3D) + 1):
            val = lautocorr3D[i - 1]
            if math.isnan(val) == True:
                dout["AUTOCORR3D" + str(i)] = "NA"
            else:
                dout["AUTOCORR3D" + str(i)] = round(val, 6)
        return dout
    else:
        for desc in _autocorr3D.keys():
            dout[desc] = "NA"
        return dout



_autocorr3D = {}
ldesc = ["AUTOCORR3D" + str(i) for i in range(1, 81)]
for desc in ldesc:
    _autocorr3D[desc] = getAUTOCORR


def GetAUTOCORR3D(mol3D):
    return getAUTOCORR(mol3D)
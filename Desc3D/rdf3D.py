from rdkit.Chem import rdMolDescriptors
import math

def getRDF(mol3D):
    dout = {}
    try: lrdf = rdMolDescriptors.CalcRDF(mol3D)
    except: lrdf = []

    if lrdf != []:
        for i in range(1, len(lrdf) + 1):
            val = lrdf[i - 1]
            if math.isnan(val) == True:
                dout["RDF" + str(i)] = "NA"
            else:
                dout["RDF" + str(i)] = round(val, 6)
        return dout

    else:
        for desc in _rdf3D.keys():
            dout[desc] = "NA"
        return dout


_rdf3D = {}
ldesc = ["RDF" + str(i) for i in range(1, 211)]
for desc in ldesc:
    _rdf3D[desc] = getRDF


def GetRDF3D(mol3D):
    return getRDF(mol3D)
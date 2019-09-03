from rdkit.Chem import rdMolDescriptors

def getRDF(mol3D):
    dout = {}
    lrdf = rdMolDescriptors.CalcRDF(mol3D)
    for i in range(1, len(lrdf) + 1):
        dout["RDF" + str(i)] = round(lrdf[i - 1], 6)
    return dout


_rdf3D = {}
ldesc = ["RDF" + str(i) for i in range(1, 211)]
for desc in ldesc:
    _rdf3D[desc] = getRDF


def GetRDF3D(mol3D):
    return getRDF(mol3D)
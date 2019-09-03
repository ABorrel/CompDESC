from rdkit.Chem import Descriptors3D, rdMolDescriptors



def rdkit3D(mol3D):
    dout = {}
    dout["Asphericity"]= Descriptors3D.Asphericity(mol3D)
    dout['Eccentricity'] = Descriptors3D.Eccentricity(mol3D)
    dout['InertialShapeFactor'] = Descriptors3D.InertialShapeFactor(mol3D)
    dout['NPR1'] = Descriptors3D.NPR1(mol3D)
    dout['NPR2'] = Descriptors3D.NPR2(mol3D)
    dout['PMI1'] = Descriptors3D.PMI1(mol3D)
    dout['PMI2'] = Descriptors3D.PMI2(mol3D)
    dout['PMI3'] = Descriptors3D.PMI3(mol3D)
    dout['RadiusOfGyration'] = Descriptors3D.RadiusOfGyration(mol3D)
    dout['SpherocityIndex'] = Descriptors3D.SpherocityIndex(mol3D)
    lautocor3D = rdMolDescriptors.CalcAUTOCORR3D(mol3D)
    for i in range(1, len(lautocor3D)+1):
        dout["AUTOCORR3D" + str(i)] = lautocor3D[i-1]
    lrdf = rdMolDescriptors.CalcRDF(mol3D)
    for i in range(1, len(lrdf)+1):
        dout["RDF" + str(i)] = lrdf[i-1]
    lmorse = rdMolDescriptors.CalcMORSE(mol3D)
    for i in range(1, len(lmorse) + 1):
        dout["MORSE" + str(i)] = lmorse[i - 1]
    lwhim = rdMolDescriptors.CalcWHIM(mol3D)
    for i in range(1, len(lwhim) + 1):
        dout["WHIM" + str(i)] = lwhim[i - 1]
    lgetaway = rdMolDescriptors.CalcGETAWAY(mol3D)
    for i in range(1, len(lgetaway) + 1):
        dout["GETAWAY" + str(i)] = lgetaway[i - 1]

    return dout
    print(dout)
import Chemical
prout = "./../trash/"

#descriptor.getLdesc("1D2D")

#SMILES = "OC(=NC1=NC2=C(N=CN2C2CC(OC(=O)C3=CC=CC=C3)C(COC(=O)C3=CC=CC=C3)O2)C(O)=N1)C1=CC=CC=C1"
#SMILES = "OO"
#SMILES = "[2H]C([2H])(C#C)N([11CH3])[C@H](C)Cc1ccccc1"
#SMILES = "NC(N)c1ccc2c(c1)[nH]c1[n+]2[Zn][n+]2c([nH]c3cc(C(N)N)ccc32)C1=O"
#SMILES = "C1=c2ccc3n2[Fe]24n5c1ccc5C=c1ccc(n12)=Cc1ccc(n14)C=3"
SMILES = "N=C(O)[C@@H](N)CS"
SMILES = "COC(C[Hg+])CNC(=O)c1ccccc1OCC(=O)O"
SMILES = "CSCC[C@H](NC(=O)[C@H](CO)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@@H](N)CO)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@@H](CCCCN)C(=O)N[C@H](C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N1CCC[C@H]1C(=O)O)C(C)C)C(C)C)C(C)C"
SMILES = "O=S(=O)(O)c1ccc([Hg])cc1"
#SMILES = "N=C(N)NCCC[C@H](N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)O"
#SMILES = "CC(C)(C)CC(C)(C)c1ccc(OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO)cc1"
SMILES = "CC(=O)N[C@H](CSSC[C@H](N)C(=O)O)C(=O)N[C@H](C)C(=O)N[C@H](CCCNC(=N)N)C(=O)N[C@H](CCCNC(=N)N)C(=O)N[C@H](CCCNC(=N)N)C(=O)N[C@H](C)C(=O)N[C@H](CCCNC(=N)N)C(N)=O"
SMILES = "CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CCCCN)C(=O)N[C@H](C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CC(C)C)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCCNC(=N)N)C(=O)O)C(C)C"

cDesc = Chemical.Chemical(SMILES, prout)
cDesc.update = 1
cDesc.prepChem()
#cDesc.computeAll2D()
#for desc in cDesc.all2D.keys():
#    print(desc, cDesc.all2D[desc])
cDesc.set3DChemical()
cDesc.computeAll3D()
for desc in cDesc.all3D.keys():
    print(desc, cDesc.all3D[desc])
sss

cDesc.prepChem()
cDesc.computePADEL2DandFP()
cDesc.computeOperaDesc()
ddd


cDesc.set3DChemical()
cDesc.computeAll3D()


pfilout = prout + "3D.csv"
filout = open(pfilout, "w")
for desc in cDesc.all3D.keys():
    filout.write("%s\t%s\n"%(desc, cDesc.all3D[desc]))
filout.close()
sss


# rdkit 3D
drdkit = rdkit3D.rdkit3D(cDesc.mol3D)
pfilout = prout + "descrdkit3D.csv"
filout = open(pfilout, "w")
for desc in drdkit.keys():
    filout.write("%s\t%s\n"%(desc, drdkit[desc]))
filout.close()
ddd





cDesc.computeAll2D()
print(len(list(cDesc.all2D.keys())))


pfilout = prout + "all.csv"
filout = open(pfilout, "w")
for desc in cDesc.all2D.keys():
    filout.write("%s\t%s\n"%(desc, cDesc.all2D[desc]))
filout.close()
ddd

drdkit = rdkitBase.computeDesc(cDesc.mol)
print("**")
print(drdkit["MaxEStateIndex"])

pfilout = prout + "descrdkit.csv"
filout = open(pfilout, "w")
for desc in drdkit.keys():
    filout.write("%s\t%s\n"%(desc, drdkit[desc]))
filout.close()


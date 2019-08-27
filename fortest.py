from Desc1D2D import rdkitBase
import descriptor


prout = "./../trash/"

SMILES = "OC(=NC1=NC2=C(N=CN2C2CC(OC(=O)C3=CC=CC=C3)C(COC(=O)C3=CC=CC=C3)O2)C(O)=N1)C1=CC=CC=C1"
#SMILES = "OO"
cDesc = descriptor.Descriptor(SMILES, "")
cDesc.computeAll2D()

pfilout = prout + "kappa.csv"
filout = open(pfilout, "w")
for desc in cDesc.kappa.keys():
    filout.write("%s\t%s\n"%(desc, cDesc.kappa[desc]))
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


# writ
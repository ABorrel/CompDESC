from Desc1D2D import rdkitBase
import descriptor



SMILES = "OC(=NC1=NC2=C(N=CN2C2CC(OC(=O)C3=CC=CC=C3)C(COC(=O)C3=CC=CC=C3)O2)C(O)=N1)C1=CC=CC=C1"
cDesc = descriptor.Descriptor(SMILES, "")

rdkitBase.computeDesc(cDesc.mol)
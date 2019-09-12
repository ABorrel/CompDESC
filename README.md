# Molecular Descriptors
Project built to harmonize Molecular Descriptor computation using Python3.6

- 8-20-19: Init project
- 8-21-19: update with RDKIT for python 3.6 consitution, molproperty
- 8-22-19: update with RDKIT for python 3.6 topological
- 8-27-19: update kappa, connectivity, bcut
- 8-30-19: Finish update for 1D and 2D descriptors
- 9-6-19: Update 3D descriptors and add also OPERA and PADEL descriptors as a supp
- 9-7-19: Install molvs in native on the conda and add in chemical class the SMILES prep process
- 9-10-19: Change table property to use the native RDKIT property
- 9-10-19: Update table of atomic property and fix bug cleanning process
- 9-12-19: Fix minor bug and optimize speed

# Dependancies
Development in python3.6 with
- RDKit (> 3.1): http://rdkit.org/docs/index.html
- molVS (> 1): https://molvs.readthedocs.io/en/latest/index.html

(additional some function will not work in case of no install)
- OPERA2.3_CL (https://github.com/kmansouri/OPERA/releases), fix the minor error in the install folder add a "/" at the path beginning
- OPERA will install PADEL in the same folder

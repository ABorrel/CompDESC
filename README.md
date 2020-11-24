# Molecular Descriptors
Project built to harmonize Molecular Descriptor computation using Python3.6

## installtion with pip
$pip install -i https://test.pypi.org/simple/ CompDesc

## Dependancies
Development in python3.6 with
- RDKit (> 3.1): http://rdkit.org/docs/index.html
- molVS (> 1): https://molvs.readthedocs.io/en/latest/index.html
- Open Babel 3.0.0 (March 2020): (http://openbabel.org/wiki/Main_Page) (sudo apt install openbabel) 
(additional some function will not work in case of no install)
- OPERA2.3_CL (https://github.com/kmansouri/OPERA/releases), fix the minor error in the install folder add a "/" at the path beginning
- OPERA will install PADEL in the same folder

## List of updates
- 8-20-19: Init project
- 8-21-19: Update with RDKIT for python 3.6 consitution, molproperty
- 8-22-19: Update with RDKIT for python 3.6 topological
- 8-27-19: Update kappa, connectivity, bcut
- 8-30-19: Finish update for 1D and 2D descriptors
- 9-6-19: Update 3D descriptors and add also OPERA and PADEL descriptors as a supp
- 9-7-19: Install molvs in native on the conda and add in chemical class the SMILES prep process
- 9-10-19: Change table property to use the native RDKIT property
- 9-10-19: Update table of atomic property and fix bug cleanning process
- 9-12-19: Fix minor bug and optimize speed
- 9-17-19: Change import and PNG file to server version
- 11-6-19: Change the path of the xml for the fp in OPERA run
- 6-9-20: Change the get list descriptors functions with OPERA desc
- 10-6-20: Add options to add salt.txt in the class and precise OS in case of molconvert
- 24-9-20: _______________ Create the package CompDesc___________________________
- 24-11-20: Add At properties in atom properties


## to do list
- check if the function getLdesc do not distrub project 
- add operating system test in the class
- add OPERA tox prediction => for now only physico-chem properties
- URGENT: FIX OPERA ERROR WITH CDK


## Usefull command lines
$python -m unittest tests/testChemical.py #unit test on Chemical class
$python setup.py sdist bdist_wheel
$python -m twine upload --repository testpypi dist/* #upload on testpypi and precise the version
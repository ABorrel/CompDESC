# List of updates
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
- 11-01-21: Minors bug and add Molecular Quantum Numbers (MQN)
- 11-01-21: Add composition descriptors
- 13-04-21: add option to build png without bg
- 20-05-21: Add fingerprint
- 24-06-21: update OPERA 2.7 and test on python 3.9
- 24-06-21: fix error with inchikey not properly generated 
- 21-11-22: Add lipinski fail -> define version 1.0
- 21-11-22: remove molvs dependency and prep SMILES with rdkit directly
- 21-11-22: add element to prep mixture
- 21-11-22: create version 1.0

# to do list
- check if the function getLdesc do not distrub project 
- add project in https://pypi.org/
- ~~remove molvs dependancy to use native rdkit with more cleaning option~~
- ~~add lipinski fail descriptor~~
- ~~add operating system test in the class~~
- ~~add OPERA tox prediction => for now only physico-chem properties~~
- ~~URGENT: FIX OPERA ERROR WITH CDK~~


# Usefull command lines
> $python -m unittest tests/testChemical.py #unit test on Chemical class

> $python setup.py sdist bdist_wheel

> $python -m twine upload --repository testpypi dist/* #upload on testpypi and precise the version
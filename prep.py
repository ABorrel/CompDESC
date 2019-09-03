from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem.SaltRemover import SaltRemover
from re import search

import toolbox

from .molVS import standardize_smiles, Standardizer

LSALT="[Co]"


LSMILESREMOVE=["[C-]#N", "[Al+3]", "[Gd+3]", "[Pt+2]", "[Au+3]", "[Bi+3]", "[Al]", "[Si+4]", "[Fe]", "[Zn]", "[Fe+2]",
               "[Ru+8]", "[Fe+]", "[Sr++]", "[Fe+3]", "[O--]", "[OH-]", "[Mn++]", "[La+3]", "[Lu+3]", "[SH-]", "[Pt+4]",
               "[Fe++]", "[W]", "[Cu+2]", "[Cr+3]", "[Tc+7]", "[Xe]", "[Tl+]", "[Zn+2]", "[F-]", "[C]", "[He]", "N#N",
               "O=O", "Cl[Ra]Cl", "[Mn+2]", "N#[N+][O-]", "II", "[Ga+3]", "[Mo+10]", "[Zn]", "[Fe]", "[Si+4]", "[Al]",
               "[B+3]"]

class prep:

    def __init__(self, chemIn, prout):
        self.chemIn = chemIn
        self.prout = prout
        self.log = ""
        self.err = 0


    def clean(self, SMICLEAN=""):

        if SMICLEAN != "":
            self.smiclean = SMICLEAN
            self.molClean = Chem.MolFromSmiles(SMICLEAN)
            return

        # Check if it is a CAS
        if not search(r"([a-z]|[A-Z])", self.chemIn):
            SMIin = self.CASDTXIDtoSMILES()
            if SMIin == "ERROR":
                self.err = 1
                return

        # Check if DSSTOX
        elif search(r"DTXSID", self.chemIn.upper()):
            SMIin = self.CASDTXIDtoSMILES()
            if SMIin == "ERROR":
                self.err = 1
                return

        else:
            SMIin = self.chemIn

        mol = Chem.MolFromSmiles(SMIin)
        s = Standardizer()

        try:
            molstandardized = s.standardize(mol)
            smilestandadized = Chem.MolToSmiles(molstandardized)
        except:
            self.log = self.log + "Normalize SMILES: ERROR DURING THE PROCESS\n"
            self.err = 1
            return


        # remove salt
        # 1.default
        remover = SaltRemover()
        mol = Chem.MolFromSmiles(smilestandadized)
        molcleandefault = remover(mol)
        # 2. Personal remover
        homeremover = SaltRemover(defnData=LSALT)
        molclean = homeremover(molcleandefault)
        smilesclean = Chem.MolToSmiles(molclean)
        # 3. SMILES remove other manual salts + fragments -> for fragment take one if exactly same compound
        lelem = smilesclean.split(".")
        if len(lelem) > 1:
            # reduce double, case of several salts are included - 255
            lelem = list(set(lelem))
            for smilesdel in LSMILESREMOVE:
                if smilesdel in lelem:
                    lelem.remove(smilesdel)
            try:
                lelem.remove("")  # case of bad smile
            except:
                pass
            if len(lelem) == 1:
                smilesclean = str(lelem[0])
            else:
                # 4. Fragments
                # Case of fragment -> stock in log file, check after to control
                self.log = self.log + "Fragments after standardization: " + smilesclean + "\n"
                self.err = 1
                return

            if smilesclean == "":
                self.log = self.log + "SMILES empty after preparation\n"
                self.err = 1
                return

        self.smiclean = smilesclean
        self.molClean = Chem.MolFromSmiles(smilesclean)


    def generateInchiKey(self):

        if not "molClean" in self.__dict__:
            self.clean()

        self.inchi = Chem.inchi.MolToInchi(self.molClean)
        self.inchikey = Chem.inchi.InchiToInchiKey(self.inchi)


    def generate3D(self, prSDF3D):

        self.prSDF3D = prSDF3D
        # generation using the method of Riniker and Landrum
        molH = Chem.AddHs(self.molClean)
        err = AllChem.EmbedMolecule(molH, AllChem.ETKDG())
        if err == 1:
            print("ERROR")#Have to do a error
        else:
            wmol = Chem.MolToMolBlock(molH)

            if not "inchikey" in self.__dict__:
                self.generateInchiKey()

            pmol = self.prSDF3D + self.inchikey + ".mol"
            fmol3D = open(pmol, "w")
            fmol3D.write(wmol)
            fmol3D.close()

            psdf3D = self.prSDF3D + self.inchikey + ".sdf"
            toolbox.babelConvertMoltoSDF(pmol, psdf3D)

            self.p3Dsdf = psdf3D


    def CASDTXIDtoSMILES(self):

        with requests.Session() as s:
            resp = s.get('https://comptox.epa.gov/dashboard/dsstoxdb/results?search=' + str(self.smiIn))
            #result = str(resp)
            resp = str(resp.text)
            try:
                SMILES = resp.split("SMILES: ")[1][0:1000]
                SMILES = SMILES.split("\\nInChI")[0]
                return SMILES
            except:
                return "ERROR"



    def parseSDFfor3DdescComputation(self):
        """
        Read the coordinates and charge of each atom in molecule from .sdf file.
        """

        dchargeSDF = {7: -3.0, 6: -2.0, 5: -1.0, 0: 0.0, 3: 1.0, 2: 2.0, 1: 3.0}  # and 4 for radical

        latoms = []

        filin = open(self.p3Dsdf, 'r')
        llines = filin.readlines()
        filin.close()

        # start at line 5 classical format
        for AtBlock in llines[4:]:
            if len(AtBlock) != 70 and len(AtBlock) != 52:
                break
            else:
                if search("IND", AtBlock):
                    continue
                X = float(AtBlock[0:10].replace(" ", ""))
                Y = float(AtBlock[10:20].replace(" ", ""))
                Z = float(AtBlock[20:30].replace(" ", ""))
                elem = AtBlock[31:34].replace(" ", "")
                charge = int(AtBlock[36:39].replace(" ", ""))
                charge = dchargeSDF[charge]

                at = [X, Y, Z, elem, charge]
                latoms.append(at)

        self.lcoords = latoms





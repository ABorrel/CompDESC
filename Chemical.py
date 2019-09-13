import functionToolbox

#2D descriptors
from Desc1D2D import constitution
from Desc1D2D import molproperty
from Desc1D2D import topology
from Desc1D2D import connectivity
from Desc1D2D import kappa
from Desc1D2D import bcut
from Desc1D2D import basak
from Desc1D2D import EStateGlobal
from Desc1D2D import moreaubroto
from Desc1D2D import moran
from Desc1D2D import geary
from Desc1D2D import charge
from Desc1D2D import moe
from Desc1D2D import morgan

#3D descriptors
from Desc3D import geo3D
from Desc3D import morse3D
from Desc3D import cpsa3D
from Desc3D import rdf3D
from Desc3D import whim3D
from Desc3D import getaway3D
from Desc3D import autocorrelation3D

#clean chemical
from prepChem import prepChem

from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem

import subprocess
from os import path, remove
from shutil import move
from re import search





class Chemical:

    # mol have to be clean before
    def __init__(self, input, prdesc):
        self.input = input
        self.err = 0
        self.log = ""
        self.prdesc = prdesc
        self.update = 0


    def prepChem(self):
        smi = prepChem.prepInput(self.input)
        print(smi, type(smi))
        if search(r"Error", smi):
            self.err = 1
            self.log = self.log + "Error: no chemical prepared\n"
        else:
            self.smiIn = smi
            smiClean = prepChem.prepSMI(smi)
            if search("Error", smiClean):
                self.err = 1
                self.log = self.log + smiClean + "\n"
            else:
                self.smi = smiClean
                self.mol = Chem.MolFromSmiles(smiClean)


    def writeSMIClean(self):
        prSMI = functionToolbox.createFolder(self.prdesc + "SMI/")
        if not "smi" in self.__dict__ and self.err == 0:
            self.prepChem()

        self.generateInchiKey()
        pfilout = prSMI + str(self.inchikey) + ".smi"
        if not path.exists(pfilout):
            filout = open(pfilout, "w")
            filout.write(self.smi)
            filout.close()
        return pfilout



    def computePNG(self):
        inchi = Chem.inchi.MolToInchi(self.mol)
        inchikey = Chem.inchi.InchiToInchiKey(inchi)

        prPNG = functionToolbox.createFolder(self.prdesc + "PNG/")
        prSMILES = functionToolbox.createFolder(self.prdesc + "SMI/")
        pSMILES = prSMILES + inchikey + ".smi"
        pPNG = prPNG + inchikey + ".png"
        if path.exists(pPNG):
            return pPNG
        else:
            if not path.exists(pSMILES):
                fSMI = open(pSMILES, "w")
                fSMI.write(str(self.smi))
                fSMI.close()
            cmd = "molconvert \"png:w500,Q100,#00000000\" " + pSMILES + " -o " + pPNG
            subprocess.Popen(cmd, shell=True)

        if path.exists(pPNG):
            return pPNG
        else:
            self.err = 1
            return "ERROR: PNG Generation"


    def generateInchiKey(self):

        if not "inchikey" in self.__dict__:
            self.inchi = Chem.inchi.MolToInchi(self.mol)
            self.inchikey = Chem.inchi.InchiToInchiKey(self.inchi)

        if self.inchikey == None:
            self.err = 1
            return 0

        return self.inchikey


    def computeAll2D(self, update = 1):

        pr2D = functionToolbox.createFolder(self.prdesc + "2D/")
        if not "inchikey" in self.__dict__:
            self.generateInchiKey()
        if self.err == 1:
            return
        pdesc2D = pr2D + self.inchikey + ".csv"
        if update == 1:
            try: remove(pdesc2D)
            except: pass
        if path.exists(pdesc2D) and path.getsize(pdesc2D) > 0:
            ddesc = functionToolbox.loadMatrixToDict(pdesc2D)
            self.all2D = ddesc
            return

        self.consti = constitution.GetConstitutional(self.mol)
        self.molprop = molproperty.GetMolecularProperty(self.mol)
        self.topology = topology.GetTopology(self.mol)
        self.connectivity = connectivity.GetConnectivity(self.mol)
        self.kappa = kappa.GetKappa(self.mol)
        self.bcut = bcut.GetBcut(self.mol)
        self.basak = basak.Getbasak(self.mol)
        self.estate = EStateGlobal.GetEState(self.mol)
        self.moreauBurto = moreaubroto.GetMBA(self.mol)
        self.moran = moran.GetMATS(self.mol)
        self.geary = geary.GetGATS(self.mol)
        self.charge = charge.GetCharge(self.mol)
        self.moe = moe.GetMOE(self.mol)
        self.morgan = morgan.GetMorgan(self.mol)

        # combine 2D
        self.all2D = {}
        self.all2D.update(deepcopy(self.consti))
        self.all2D.update(deepcopy(self.molprop))
        self.all2D.update(deepcopy(self.topology))
        self.all2D.update(deepcopy(self.connectivity))
        self.all2D.update(deepcopy(self.kappa))
        self.all2D.update(deepcopy(self.bcut))
        self.all2D.update(deepcopy(self.basak))
        self.all2D.update(deepcopy(self.estate))
        self.all2D.update(deepcopy(self.moreauBurto))
        self.all2D.update(deepcopy(self.moran))
        self.all2D.update(deepcopy(self.geary))
        self.all2D.update(deepcopy(self.charge))
        self.all2D.update(deepcopy(self.moe))
        self.all2D.update(deepcopy(self.morgan))



    def set3DChemical(self, psdf3D = ""):

        if "psdf3D" in self.__dict__ or psdf3D != "":
            self.coords = functionToolbox.parseSDFfor3DdescComputation(psdf3D)
            return

        else:
            prSDF3D = functionToolbox.createFolder(self.prdesc + "SDF3D/")
            prMOLCLEAN = functionToolbox.createFolder(self.prdesc + "MOLCLEAN/")

            # control if exist
            if not "inchikey" in self.__dict__:
                self.generateInchiKey()
            pmol = prMOLCLEAN + self.inchikey + ".mol"
            psdf3D = prSDF3D + self.inchikey + ".sdf"

            if self.update == 1:
                try: remove(pmol)
                except: pass
                try:remove(psdf3D)
                except: pass

            if path.exists(pmol) and path.exists(psdf3D):
                self.coords = functionToolbox.parseSDFfor3DdescComputation(psdf3D)
                self.psdf3D = psdf3D
                self.mol3D = Chem.MolFromMolFile(pmol)
                return

            # have to generate the 3D
            molH = Chem.AddHs(self.mol)
            err = AllChem.EmbedMolecule(molH, AllChem.ETKDG())
            #err = AllChem.MMFFOptimizeMolecule(molH)
            if err == 1 :#or err == -1:
                self.err = 1
                print("ERROR 3D generation")  # Have to do a error
                return

            self.mol3D = molH

            # to write
            wmol = Chem.MolToMolBlock(molH)
            fmol3D = open(pmol, "w")
            fmol3D.write(wmol)
            fmol3D.close()

            functionToolbox.babelConvertMoltoSDF(pmol, psdf3D, window=0, update=self.update)

            self.coords = functionToolbox.parseSDFfor3DdescComputation(psdf3D)
            self.psdf3D = psdf3D


    def computeAll3D(self, update = 1):
        self.update = update
        pr3D = functionToolbox.createFolder(self.prdesc + "3D/")

        # control if error of 3D generation
        if self.err == 1:
            return

        if not "inchikey" in self.__dict__:
            self.generateInchiKey()
        pdesc3D = pr3D + self.inchikey + ".csv"
        if self.update == 1:
            try: remove(pdesc3D)
            except: pass

        if path.exists(pdesc3D) and path.getsize(pdesc3D) > 0:
            ddesc = functionToolbox.loadMatrixToDict(pdesc3D)
            self.all3D = ddesc
            return

        # compute descriptors
        self.geo3D = geo3D.GetGeo3D(self.coords, self.mol3D)
        self.morse3D = morse3D.GetMorse3D(self.mol3D)
        self.whim3D = whim3D.GetWHIM3D(self.mol3D)
        self.rdf3D = rdf3D.GetRDF3D(self.mol3D)
        self.getaway3D = getaway3D.GetGETAWAY3D(self.mol3D)
        self.autocorr3D = autocorrelation3D.GetAUTOCORR3D(self.mol3D)
        self.CPSA3D = cpsa3D.GetCPSA3D(self.coords)

        # combine 3D
        self.all3D = {}
        self.all3D.update(deepcopy(self.geo3D))
        self.all3D.update(deepcopy(self.CPSA3D))
        self.all3D.update(deepcopy(self.rdf3D))
        self.all3D.update(deepcopy(self.morse3D))
        self.all3D.update(deepcopy(self.whim3D))
        self.all3D.update(deepcopy(self.getaway3D))
        self.all3D.update(deepcopy(self.autocorr3D))



    def computePADEL2DandFP(self, PPADEL=""):
        prPadelDesc = functionToolbox.createFolder(self.prdesc + "PADEL_desc/")
        prPadelFp = functionToolbox.createFolder(self.prdesc + "PADEL_fp/")
        prPadelTemp = functionToolbox.createFolder(self.prdesc + "PADEL_temp/", 1)
        if "smi" in self.__dict__:
            if not "inchikey" in self.__dict__:
                self.generateInchiKey()
            ppadel_desc = prPadelDesc + self.inchikey + ".csv"
            ppadel_FP = prPadelFp + self.inchikey + ".csv"
            if path.exists(ppadel_desc) and path.exists(ppadel_FP):
                self.ppadel_desc = ppadel_desc
                self.ppadel_FP = ppadel_FP
            else:
                pSMI = prPadelTemp + self.inchikey + ".smi"
                fSMI = open(pSMI, "w")
                fSMI.write(self.smi)
                fSMI.close()
                pdesc = functionToolbox.runPadelDesc(prPadelTemp, PPADEL)
                move(pdesc, ppadel_desc)
                pFP = functionToolbox.runPadelFP(prPadelTemp, PPADEL)
                move(pFP, ppadel_FP)
                self.ppadel_desc = ppadel_desc
                self.ppadel_FP = ppadel_FP


    def computeOperaDesc(self, POPERA = "", PMATLAB = ""):

        if not "ppadel_desc" in self.__dict__:
            self.computePADEL2DandFP()

        prOPERA = functionToolbox.createFolder(self.prdesc + "OPERA/")
        pfilout = prOPERA + self.inchikey + ".csv"
        functionToolbox.runOPERA(self.ppadel_desc, self.ppadel_FP, pfilout, POPERA, PMATLAB)
        self.pOPERA = pfilout


    def writeSDF(self, psdf, name):

        if not "mol" in self.__dict__:
            print("No mol load")
            return 1

        if path.exists(psdf) and path.getsize(psdf) > 100:
            fsdf = open(psdf, "a")
            fsdf.write("\n$$$$\n")
        else:
            fsdf = open(psdf, "w")

        molH = Chem.AddHs(self.mol)
        fsdf.write(str(name) + "\n" + Chem.MolToMolBlock(molH)[1:])
        fsdf.close()


    def writeMatrix(self, typedesc):
        if typedesc == "2D":
            pr2D = functionToolbox.createFolder(self.prdesc + "2D/")
            if "all2D" in self.__dict__:
                filin = open(pr2D + self.inchikey + ".csv", "w")
                filin.write("%s\n"%("\t".join(self.all2D.keys())))
                filin.write("%s\n"%("\t".join([str(self.all2D[k]) for k in self.all2D.keys()])))
                filin.close()

        if typedesc == "3D":
            pr3D = functionToolbox.createFolder(self.prdesc + "3D/")
            if "all3D" in self.__dict__:
                filin = open(pr3D + self.inchikey + ".csv", "w")
                filin.write("%s\n"%("\t".join(self.all3D.keys())))
                filin.write("%s\n"%("\t".join([str(self.all3D[k]) for k in self.all3D.keys()])))
                filin.close()



### Load a list of descriptors ###
##################################
def getLdesc (typeDesc):

    lout = []
    if typeDesc == "1D2D" or typeDesc == "1D" or typeDesc == "2D" or typeDesc == "all":
        lout = lout + list(constitution._constitutional.keys()) + list(molproperty._molProperty.keys()) + list(topology._topology.keys()) +\
                list(connectivity._connectivity.keys()) + list(kappa._kappa.keys()) + list(bcut._bcut.keys()) + list(basak._basak.keys()) +\
                list(EStateGlobal._EState.keys()) + list(moreaubroto._MBA.keys()) + list(moran._moran.keys()) + list(geary._geary.keys()) +\
                list(charge._charge.keys()) + list(moe._moe.keys()) + list(morgan._morgan.keys())

    elif typeDesc == "3D" or typeDesc == "all":
        lout = lout + list(autocorrelation3D._autocorr3D.keys()) + list(cpsa3D._cpsa3D.keys()) + list(geo3D._geo3D.keys()) + \
            list(getaway3D._getaway3D.keys()) + list(morse3D._morse3D.keys()) + list(rdf3D._rdf3D.keys()) + list(whim3D._whim3D.keys())

    return lout
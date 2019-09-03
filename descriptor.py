import toolbox

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

#from rdkit.Chem import Descriptors
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdPartialCharges as GMCharge

import subprocess
from os import path, remove


def getLdesc (typeDesc):

    lout = []
    if typeDesc == "1D2D" or typeDesc == "all":
        lout = lout + list(constitution._constitutional.keys()) + list(molproperty._molProperty.keys()) + list(topology._topology.keys()) +\
                list(connectivity._connectivity.keys()) + list(kappa._kappa.keys()) + list(bcut._bcut.keys()) + list(basak._basak.keys()) +\
                list(EStateGlobal._EState.keys()) + list(moreaubroto._MBA.keys()) + list(moran._moran.keys()) + list(geary._geary.keys()) +\
                list(charge._charge.keys()) + list(moe._moe.keys()) + list(morgan._morgan.keys())

    elif typeDesc == "3D" or typeDesc == "all":
        lout = lout + list(autocorrelation3D._autocorr3D.keys()) + list(cpsa3D._cpsa3D.keys()) + list(geo3D._geo3D.keys()) + \
            list(getaway3D._getaway3D.keys()) + list(morse3D._morse3D.keys()) + list(rdf3D._rdf3D.keys()) + list(whim3D._whim3D.keys())

    return lout



class Descriptor:

    # mol have to be clean before
    def __init__(self, SMICLEAN, prdesc):
        self.smi = SMICLEAN
        self.mol = Chem.MolFromSmiles(SMICLEAN)
        self.err = 0
        self.prdesc = prdesc


    def computePNG(self, prSMILES, prPNG):

        inchi = Chem.inchi.MolToInchi(self.mol)
        inchikey = Chem.inchi.InchiToInchiKey(inchi)

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

        self.inchi = Chem.inchi.MolToInchi(self.mol)
        self.inchikey = Chem.inchi.InchiToInchiKey(self.inchi)


    def computeAll2D(self, update = 1):
        if path.exists(self.prdesc + "desc2D.txt") and update == 0:
            if path.getsize(self.prdesc + "desc2D.txt") > 100:
                ddesc = toolbox.loadMatrixToDict(self.prdesc + "desc2D.txt")
                self.all2D = ddesc
                return
            else:
                remove(self.prdesc + "_2D.txt")

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
        if psdf3D == "":
            prSDF3D = toolbox.createFolder(self.prdesc + "SDF3D/")
            prMOLCLEAN = toolbox.createFolder(self.prdesc + "MOLCLEAN/")
            # have to generate the 3D
            molH = Chem.AddHs(self.mol)
            err = AllChem.EmbedMolecule(molH, AllChem.ETKDG())
            print(err)
            if err == 1:
                self.err = 1
                print("ERROR 3D generation")  # Have to do a error
                return

            self.mol3D = molH

            # to write
            wmol = Chem.MolToMolBlock(molH)
            if not "inchikey" in self.__dict__:
                self.generateInchiKey()

            pmol = prMOLCLEAN + self.inchikey + ".mol"
            fmol3D = open(pmol, "w")
            fmol3D.write(wmol)
            fmol3D.close()

            psdf3D = prSDF3D + self.inchikey + ".sdf"
            toolbox.babelConvertMoltoSDF(pmol, psdf3D)

        self.coords = toolbox.parseSDFfor3DdescComputation(psdf3D)


    def computeAll3D(self):

        if path.exists(self.prdesc + "_3D.txt"):
            if path.getsize(self.prdesc + "_3D.txt") > 100:
                self.all3D = toolbox.loadMatrixToDict(self.prdesc + "_3D.txt")
                return
            else:
                remove(self.prdesc + "_3D.txt")

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



    def writeMatrix(self, typedesc):
        if typedesc == "2D":
            if "all2D" in self.__dict__:
                filin = open(self.prdesc + "_2D.txt", "w")
                filin.write("%s\n"%("\t".join(self.all2D.keys())))
                filin.write("%s\n"%("\t".join([str(self.all2D[k]) for k in self.all2D.keys()])))
                filin.close()

        if typedesc == "3D":
            if "all3D" in self.__dict__:
                filin = open(self.prdesc + "_3D.txt", "w")
                filin.write("%s\n"%("\t".join(self.all3D.keys())))
                filin.write("%s\n"%("\t".join([str(self.all3D[k]) for k in self.all3D.keys()])))
                filin.close()



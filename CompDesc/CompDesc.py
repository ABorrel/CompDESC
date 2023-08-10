import sys
from . import functionToolbox

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
from Desc1D2D import MQNs

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
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit import DataStructs

from os import path, remove, system, name
from shutil import move
from re import search
from random import randint
from shutil import rmtree

class CompDesc:

    # mol have to be clean before
    def __init__(self, input, prdesc, update=0, p_salts = ""):
        self.input = input
        self.err = 0
        self.log = ""
        self.prdesc = prdesc
        self.update = 0
        self.p_salts = p_salts
        self.isfrag = 0 #export mixture or frag
        #p_repertory = str(pathlib.Path(__file__).parent.absolute())
        #self.p_xml = p_repertory + "/desc_fp.xml"

    def prepChem(self):
        smi = prepChem.prepInput(self.input)
        if search(r"Error", smi):
            self.err = 1
            self.log = self.log + "Error: no chemical prepared\n"
        else:
            self.smiIn = smi
            smiClean = prepChem.prepSMI(smi)
            if search("Error", smiClean):
                self.err = 1
                self.log = self.log + smiClean + "\n"
            elif search("Mixture", smiClean):
                self.err = 1
                self.isfrag = 1
                self.smi = smiClean.split(": ")[-1] # keep in string to not have to make too much change
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

    def computePNG(self, prPNG = "", bg="white"):
        if self.err == 1:
            return "ERROR: PNG Generation"

        if not "inchikey" in self.__dict__:
            self.generateInchiKey()
        if self.err == 1:
            return "ERROR: PNG Generation"

        if prPNG == "":
            prPNG = functionToolbox.createFolder(self.prdesc + "PNG/")
            prSMILES = functionToolbox.createFolder(self.prdesc + "SMI/")
        else:
            prSMILES = functionToolbox.createFolder(prPNG + "SMI/")

        pSMILES = prSMILES + self.inchikey + ".smi"
        pPNG = prPNG + self.inchikey + ".png"
        if path.exists(pPNG):
            return pPNG
        else:
            if not path.exists(pSMILES):
                fSMI = open(pSMILES, "w")
                fSMI.write(str(self.smi))
                fSMI.close()
            if name == "nt":
                if bg == "white":
                    cmd = "C:/\"Program Files\"/ChemAxon/MarvinSuite/bin/molconvert \"png:w500,Q100\" " + pSMILES + " -o " + pPNG
                else:
                    cmd = "C:/\"Program Files\"/ChemAxon/MarvinSuite/bin/molconvert \"png:w500,Q100,#00000000\" " + pSMILES + " -o " + pPNG
            else:
                if bg == "white":
                    cmd = "molconvert \"png:w500,Q100\" " + pSMILES + " -o " + pPNG
                else:
                    cmd = "molconvert \"png:w500,Q100,#00000000\" " + pSMILES + " -o " + pPNG
            system(cmd)
            #subprocess.Popen(cmd, shell=True)

        if path.exists(pPNG):
            return pPNG
        else:
            self.err = 1
            return "ERROR: PNG Generation"

    def generateInchiKey(self):

        if not "inchikey" in self.__dict__:
            if not "mol" in self.__dict__:
                self.prepChem()
            if self.err == 0:
                try:
                    self.inchi = Chem.inchi.MolToInchi(self.mol)
                    self.inchikey = Chem.inchi.InchiToInchiKey(self.inchi)
                except:
                    self.err = 1
                    return 0
            else:
                self.err = 1
                return 0

        if self.inchikey == None:
            self.err = 1
            return 0

        return self.inchikey

    def computeAll2D(self):

        if not "inchikey" in self.__dict__:
            self.generateInchiKey()
        
        if self.err == 1: # from inchkey
            return

        pr2D = self.prdesc + "2D/"
        pdesc2D = pr2D + self.inchikey + ".csv"

        if self.update == 1:
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
        self.mqns = MQNs.GetMQNs(self.mol)

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
        self.all2D.update(deepcopy(self.mqns))

    def set3DChemical(self, psdf3D = "", w=0):

        if "psdf3D" in self.__dict__ or psdf3D != "":
            self.coords = functionToolbox.parseSDFfor3DdescComputation(psdf3D)
            return

        else:
            prSDF3D = functionToolbox.createFolder(self.prdesc + "SDF3D/")
            prMOLCLEAN = functionToolbox.createFolder(self.prdesc + "MOLCLEAN/")

            # control if exist
            if not "inchikey" in self.__dict__:
                self.generateInchiKey()
            if self.err == 1:
                return
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
            err = AllChem.EmbedMolecule(molH,randomSeed=0xf00d)
            if err == 1 :#or err == -1:
                self.err = 1
                print("ERROR 3D generation")  # Have to do a error
                return

            Chem.MolToMolBlock(molH)
            self.mol3D = molH

            # write native rdkit
            if w == 1:
                with Chem.SDWriter(psdf3D) as w:
                    w.write(self.mol3D)

            self.coords = functionToolbox.parseSDFfor3DdescComputation(Chem.MolToMolBlock(self.mol3D))
            # add second control of 3d computation
            l_y = [atom[2] for atom in self.coords]
            if l_y[0] == 0.0 and l_y[-1] == 0.0:
                self.err = 1 # case pb in the 3D
            else:
                self.psdf3D = psdf3D

    def computeAll3D(self):
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

    def computePADEL2DFPandCDK(self, PPADEL="", PCDK=""):
        prPadelDesc = functionToolbox.createFolder(self.prdesc + "PADEL_desc/")
        prPadelFp = functionToolbox.createFolder(self.prdesc + "PADEL_fp/")
        prCDK = functionToolbox.createFolder(self.prdesc + "cdk_desc/")
        a = randint(0,100000)
        prPadelTemp = functionToolbox.createFolder(self.prdesc + "PADEL_" + str(a) + "/", 1)
        prCDKTemp = functionToolbox.createFolder(self.prdesc + "CDK_" + str(a) + "/", 1)
        if "smi" in self.__dict__:
            if not "inchikey" in self.__dict__:
                self.generateInchiKey()
            ppadel_desc = prPadelDesc + self.inchikey + ".csv"
            ppadel_FP = prPadelFp + self.inchikey + ".csv"
            pcdk_desc = prCDK + self.inchikey + ".csv"
            if path.exists(ppadel_desc) and path.exists(ppadel_FP):# and path.exists(pcdk_desc):
                self.ppadel_desc = ppadel_desc
                self.ppadel_FP = ppadel_FP
                self.pcdk_desc = pcdk_desc
            else:
                # for PADEL
                pSMI = prPadelTemp + self.inchikey + ".smi"
                fSMI = open(pSMI, "w")
                fSMI.write(self.smi)
                fSMI.close()

                # for CDK
                pSMI_CDK = prCDKTemp + self.inchikey + ".smi"
                fSMI_CDK = open(pSMI_CDK, "w")
                fSMI_CDK.write(self.smi)
                fSMI_CDK.close()

                #desc PADEL
                pdesc = functionToolbox.runPadelDesc(prPadelTemp, PPADEL)
                move(pdesc, ppadel_desc)
                #FP PADEL 
                pFP = functionToolbox.runPadelFP(prPadelTemp, PPADEL)#, self.p_xml)
                move(pFP, ppadel_FP)
                #desc CDK
                pdescCDK =  functionToolbox.runCDKDesc(pSMI_CDK, prCDKTemp, PCDK)
                move(pdescCDK, pcdk_desc)

                self.ppadel_desc = ppadel_desc
                self.ppadel_FP = ppadel_FP
                self.pcdk_desc = pcdk_desc
        try:         
            try: 
                remove(prPadelTemp)
                remove(prCDKTemp)
            except: 
                rmtree(prPadelTemp)
                rmtree(prCDKTemp)
        except: pass

    def computeOperaDesc(self, POPERA = "", PMATLAB = "", onlyPhysChem=0):
        """
        Function used on the server => need to check cdk
        """

        if not "ppadel_desc" in self.__dict__:
            self.computePADEL2DFPandCDK()

        prOPERA = functionToolbox.createFolder(self.prdesc + "OPERA/")
        pfilout = prOPERA + self.inchikey + ".csv"

        #functionToolbox.runOPERA(self.ppadel_desc, self.ppadel_FP, pfilout, POPERA, PMATLAB)
        functionToolbox.runOPERA(self.ppadel_desc, self.ppadel_FP, self.pcdk_desc, pfilout, POPERA, PMATLAB, onlyPhysChem)
        self.pOPERA = pfilout
        if path.exists(self.pOPERA) and path.getsize(self.pOPERA) > 0:
            d_opera = functionToolbox.loadMatrixToDict(self.pOPERA, ",")
            self.allOPERA = d_opera
        else:
            self.allOPERA = {}
            self.err = 1
        
    def computeOPERAFromChem(self, POPERA = "", PMATLAB = "", update=0):

        if not "inchikey" in self.__dict__:
                self.generateInchiKey()

        prOPERA = functionToolbox.createFolder(self.prdesc + "OPERA/")
        pfilout = prOPERA + self.inchikey + ".csv"

        psmi = prOPERA + self.inchikey + ".smi"
        fsmi = open(psmi, "w")
        fsmi.write(self.smi)
        fsmi.close()

        functionToolbox.runOPERAFromChem(psmi, pfilout, POPERA, PMATLAB, update=update)
        self.pOPERA = pfilout

        d_opera = functionToolbox.loadMatrixToDict(self.pOPERA, ",")
        self.allOPERA = d_opera

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

    def getLdesc (self, typeDesc):
        """Get list of descriptor"""
        lout = []
        if typeDesc == "1D2D" or typeDesc == "1D" or typeDesc == "2D" or typeDesc == "all":
            lout = lout + list(constitution._constitutional.keys()) + list(molproperty._molProperty.keys()) + list(topology._topology.keys()) +\
                    list(connectivity._connectivity.keys()) + list(kappa._kappa.keys()) + list(bcut._bcut.keys()) + list(basak._basak.keys()) +\
                    list(EStateGlobal._EState.keys()) + list(moreaubroto._MBA.keys()) + list(moran._moran.keys()) + list(geary._geary.keys()) +\
                    list(charge._charge.keys()) + list(moe._moe.keys()) + list(morgan._morgan.keys()) + list(MQNs._mqn.keys())

        if typeDesc == "3D" or typeDesc == "all":
            lout = lout + list(autocorrelation3D._autocorr3D.keys()) + list(cpsa3D._cpsa3D.keys()) + list(geo3D._geo3D.keys()) + \
                list(getaway3D._getaway3D.keys()) + list(morse3D._morse3D.keys()) + list(rdf3D._rdf3D.keys()) + list(whim3D._whim3D.keys())

        if typeDesc == "OPERA" or typeDesc == "all":
            lout = lout + ["MolWeight", "nbAtoms", "nbHeavyAtoms", "nbC", "nbO", "nbN" ,"nbAromAtom","nbRing","nbHeteroRing","Sp3Sp2HybRatio","nbRotBd","nbHBdAcc","ndHBdDon","nbLipinskiFailures","TopoPolSurfAir","MolarRefract","CombDipolPolariz","LogP_pred","MP_pred","BP_pred","LogVP_pred","LogWS_pred","LogHL_pred","RT_pred","LogKOA_pred","ionization","pKa_a_pred","pKa_b_pred","LogD55_pred","LogD74_pred","LogOH_pred","LogBCF_pred","BioDeg_LogHalfLife_pred","ReadyBiodeg_pred","LogKM_pred","LogKoc_pred", "FUB_pred", "Clint_pred"]
        
        if typeDesc == "KNIME":
            lout = lout + ["SlogP","SMR","LabuteASA","TPSA","AMW","ExactMW","NumLipinskiHBA","NumLipinskiHBD","NumRotatableBonds","NumHBD","NumHBA","NumAmideBonds","NumHeteroAtoms","NumHeavyAtoms","NumAtoms","NumStereocenters","NumUnspecifiedStereocenters","NumRings","NumAromaticRings","NumSaturatedRings","NumAliphaticRings","NumAromaticHeterocycles","NumSaturatedHeterocycles","NumAliphaticHeterocycles","NumAromaticCarbocycles","NumSaturatedCarbocycles","NumAliphaticCarbocycles","FractionCSP3","Chi0v","Chi1v","Chi2v","Chi3v","Chi4v","Chi1n","Chi2n","Chi3n","Chi4n","HallKierAlpha","kappa1","kappa2","kappa3","slogp_VSA1","slogp_VSA2","slogp_VSA3","slogp_VSA4","slogp_VSA5","slogp_VSA6","slogp_VSA7","slogp_VSA8","slogp_VSA9","slogp_VSA10","slogp_VSA11","slogp_VSA12","smr_VSA1","smr_VSA2","smr_VSA3","smr_VSA4","smr_VSA5","smr_VSA6","smr_VSA7","smr_VSA8","smr_VSA9","smr_VSA10","peoe_VSA1","peoe_VSA2","peoe_VSA3","peoe_VSA4","peoe_VSA5","peoe_VSA6","peoe_VSA7","peoe_VSA8","peoe_VSA9","peoe_VSA10","peoe_VSA11","peoe_VSA12","peoe_VSA13","peoe_VSA14","MQN1","MQN2","MQN3","MQN4","MQN5","MQN6","MQN7","MQN8","MQN9","MQN10","MQN11","MQN12","MQN13","MQN14","MQN15","MQN16","MQN17","MQN18","MQN19","MQN20","MQN21","MQN22","MQN23","MQN24","MQN25","MQN26","MQN27","MQN28","MQN29","MQN30","MQN31","MQN32","MQN33","MQN34","MQN35","MQN36","MQN37","MQN38","MQN39","MQN40","MQN41","MQN42"] 

        return lout

    def computeOperaFromListSMI(self, pfilout):

        functionToolbox.runOPERAFromSmi(self.input, pfilout=pfilout, update=self.update)

    def convertDesc2DtoKnimeDesc(self):

        l_desc_KNIME = ["SlogP","SMR","LabuteASA","TPSA","AMW","ExactMW","NumLipinskiHBA","NumLipinskiHBD","NumRotatableBonds","NumHBD","NumHBA","NumAmideBonds","NumHeteroAtoms","NumHeavyAtoms","NumAtoms","NumStereocenters","NumUnspecifiedStereocenters","NumRings","NumAromaticRings","NumSaturatedRings","NumAliphaticRings","NumAromaticHeterocycles","NumSaturatedHeterocycles","NumAliphaticHeterocycles","NumAromaticCarbocycles","NumSaturatedCarbocycles","NumAliphaticCarbocycles","FractionCSP3","Chi0v","Chi1v","Chi2v","Chi3v","Chi4v","Chi1n","Chi2n","Chi3n","Chi4n","HallKierAlpha","kappa1","kappa2","kappa3","slogp_VSA1","slogp_VSA2","slogp_VSA3","slogp_VSA4","slogp_VSA5","slogp_VSA6","slogp_VSA7","slogp_VSA8","slogp_VSA9","slogp_VSA10","slogp_VSA11","slogp_VSA12","smr_VSA1","smr_VSA2","smr_VSA3","smr_VSA4","smr_VSA5","smr_VSA6","smr_VSA7","smr_VSA8","smr_VSA9","smr_VSA10","peoe_VSA1","peoe_VSA2","peoe_VSA3","peoe_VSA4","peoe_VSA5","peoe_VSA6","peoe_VSA7","peoe_VSA8","peoe_VSA9","peoe_VSA10","peoe_VSA11","peoe_VSA12","peoe_VSA13","peoe_VSA14","MQN1","MQN2","MQN3","MQN4","MQN5","MQN6","MQN7","MQN8","MQN9","MQN10","MQN11","MQN12","MQN13","MQN14","MQN15","MQN16","MQN17","MQN18","MQN19","MQN20","MQN21","MQN22","MQN23","MQN24","MQN25","MQN26","MQN27","MQN28","MQN29","MQN30","MQN31","MQN32","MQN33","MQN34","MQN35","MQN36","MQN37","MQN38","MQN39","MQN40","MQN41","MQN42"] 
       
        if not "all2D" in list(self.__dict__.keys()):
            self.err = 1
            print("No descriptors computed")
            return 
        else:
            d_desc_out = {}
            i = 0
            imax = len(l_desc_KNIME)
            while i < imax:
                desc = l_desc_KNIME[i]
                if search("peoe_", desc):
                    d_desc_out[desc] = self.all2D[desc.upper()]
                elif search("slogp_", desc):
                    d_desc_out[desc] = self.all2D["SlogP_" + desc.split("_")[-1]]
                elif search("smr_", desc):
                    d_desc_out[desc] = self.all2D[desc.upper()]
                elif search("kappa", desc):
                    d_desc_out[desc] = self.all2D["s%s"%(desc)]
                elif search("Chi.v", desc):
                    d_desc_out[desc] = self.all2D["Chi%s%s"%(desc[4], desc[3])]  
                elif search("Chi.n", desc):
                    d_desc_out[desc] = self.all2D["Chiv%s"%(desc[3])] 
                elif desc == 'ExactMW':
                    d_desc_out[desc] = self.all2D["ExactMolWt"] 
                elif desc == 'NumAtoms':
                    d_desc_out[desc] = self.all2D["NumAllatoms"]
                elif desc == 'NumHeteroAtoms':
                    d_desc_out[desc] = self.all2D["NumHeteroatoms"]
                elif desc == 'NumHeavyAtoms':
                    d_desc_out[desc] = self.all2D["HeavyAtomCount"]
                elif desc == 'NumRings':
                    d_desc_out[desc] = self.all2D["RingCount"]
                elif desc == 'NumHBD':
                    d_desc_out[desc] = self.all2D["NumHDonors"]
                elif desc == 'NumHBA':
                    d_desc_out[desc] = self.all2D["NumHAcceptors"]
                elif desc == 'SlogP':
                    d_desc_out[desc] = self.all2D["MolLogP"]
                elif desc == 'SMR':
                    d_desc_out[desc] = self.all2D["MolMR"]
                elif desc == 'AMW':
                    d_desc_out[desc] = self.all2D["MolWt"]
                else:
                    d_desc_out[desc] = self.all2D[desc]
                i = i + 1
            
            # reload all2d descripor
            self.all2D = d_desc_out
                
    def computeFP(self, typeFP):

        if not "mol" in self.__dict__:
            self.log = self.log + "No smiles prepared\n"
            self.err = 1
        else:
            d_FP = {}
            if typeFP == "Mol" or typeFP == "All":
                d_FP["Mol"] = FingerprintMols.FingerprintMol(self.mol)
            if typeFP == "MACCS" or typeFP == "All":
                d_FP["MACCS"] = MACCSkeys.GenMACCSKeys(self.mol)
            if typeFP == "pairs" or typeFP == "All":
                d_FP["pairs"] = Pairs.GetAtomPairFingerprint(self.mol)
            if typeFP == "Torsion" or typeFP == "All":
                d_FP["Torsion"] = Torsions.GetTopologicalTorsionFingerprint(self.mol)
            if typeFP == "Morgan" or typeFP == "All":
                d_FP["Morgan"] = AllChem.GetMorganFingerprint(self.mol, 2)
            
            self.d_FP = d_FP

    def computeSimilarityFP(self, c_chem, typeFP, typeMetric):
        
        try:
            if typeMetric == 'Tanimoto':
                return DataStructs.TanimotoSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "Dice":
                return DataStructs.DiceSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "Cosine":
                return DataStructs.CosineSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "Sokal":
                return DataStructs.SokalSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "Russel":
                return DataStructs.RusselSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "RogotGoldberg":
                return DataStructs.RogotGoldbergSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "AllBit":
                return DataStructs.AllBitSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "Kulczynski":
                return DataStructs.KulczynskiSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "McConnaughey":
                return DataStructs.McConnaugheySimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "Asymmetric":
                return DataStructs.AsymmetricSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
            elif typeMetric == "BraunBlanquet":
                return DataStructs.BraunBlanquetSimilarity(self.d_FP[typeFP], c_chem.d_FP[typeFP])
        except:
            print("Combination %s and %s not supported"%(typeFP, typeMetric))
            self.log = "%sCombination %s and %s not supported\n"%(self.log, typeFP, typeMetric)
            return "NA"


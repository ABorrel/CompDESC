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
#from .geo3D import GetGeo3D
#from .cpsa3D import GetCPSA3D
#from .rdf3D import GetRdf3D
#from .morse3D import GetMorse3D
#from .whim3D import GetWhim3D


#from rdkit.Chem import Descriptors
from copy import deepcopy
from rdkit import Chem


import subprocess
from os import path, remove



L3D = ['RDFC6', 'MoRSEN11', 'RDFU8', 'RDFU9', 'RDFU2', 'RDFU3', 'MoRSEN5', 'RDFU1', 'RDFU6', 'RDFU7', 'RDFU4', 'RDFU5',
       'Harary3D', 'P2u', 'MoRSEM6', 'MoRSEM7', 'MoRSEM4', 'MoRSEM5', 'MoRSEM2', 'MoRSEM3', 'MoRSEE30', 'MoRSEM1',
       'MoRSEN4', 'MoRSEM8', 'MoRSEM9', 'MoRSEU10', 'MoRSEU11', 'MoRSEU12', 'MoRSEU13', 'MoRSEU14', 'MoRSEU15',
       'MoRSEU16', 'MoRSEU17', 'MoRSEU18', 'MoRSEU19', 'FPSA3', 'FPSA2', 'FPSA1', 'GeDi', 'MoRSEV19', 'MoRSEN10',
       'MoRSEV13', 'SPAN', 'MoRSEV11', 'MoRSEV10', 'MoRSEV17', 'MoRSEV16', 'MoRSEV15', 'MoRSEN16', 'RDFM14', 'RDFM15',
       'RDFM16', 'RDFE9', 'RDFM10', 'RDFM11', 'RDFM12', 'RASA', 'RDFE2', 'RDFE3', 'MoRSEC25', 'MoRSEC24', 'RDFM18',
       'RDFM19', 'grav', 'RDFE5', 'WNSA1', 'WNSA2', 'WNSA3', 'L2p', 'RDFP15', 'RDFP14', 'RDFP17', 'RDFP16', 'RDFP11',
       'RDFP10', 'RDFP13', 'RDFP12', 'MoRSEP30', 'RDFP19', 'RDFP18', 'E2p', 'Dm', 'P3e', 'MoRSEM18', 'MoRSEM19',
       'Petitj3D', 'MoRSEM10', 'MoRSEM11', 'MoRSEM12', 'MoRSEM13', 'MoRSEM14', 'MoRSEM15', 'MoRSEM16', 'MoRSEM17',
       'RDFC27', 'RDFC26', 'RDFC25', 'RDFC24', 'RDFC23', 'RDFC22', 'RDFC21', 'RDFC20', 'MoRSEC28', 'RDFU30', 'RDFC29',
       'RDFC28', 'MoRSEU30', 'L1u', 'L1v', 'L2v', 'L1p', 'RDFP5', 'RDFP4', 'RDFP7', 'RDFP6', 'RDFP1', 'RDFP3', 'RDFP2',
       'L1e', 'RDFP9', 'RDFP8', 'MoRSEP5', 'P1e', 'MoRSEP4', 'PSA', 'MoRSEP7', 'P1p', 'MoRSEP6', 'RDFE18', 'RDFE19',
       'P1v', 'RDFE14', 'RDFE15', 'RDFE16', 'RDFE17', 'PPSA1', 'RDFE11', 'RDFE12', 'PPSA2', 'MoRSEP11', 'MoRSEP10',
       'MoRSEP13', 'MoRSEP12', 'RPCS', 'MoRSEP14', 'MoRSEN9', 'MoRSEN8', 'DPSA1', 'MoRSEC30', 'DPSA3', 'DPSA2',
       'MoRSEN3', 'MoRSEN2', 'MoRSEN1', 'RDFP30', 'E2e', 'MoRSEN17', 'L3e', 'TASA', 'RDFC19', 'MoRSEV14', 'MoRSEM30',
       'MoRSEP8', 'L3v', 'RDFC16', 'L3u', 'RDFV30', 'L3p', 'RDFC14', 'W3DH', 'RDFC15', 'MoRSEC23', 'MoRSEN15',
       'MoRSEP16', 'RPSA', 'P3m', 'MEcc', 'MoRSEC22', 'MoRSEN14', 'MoRSEP1', 'MoRSEN23', 'P3p', 'P3v', 'MoRSEP19',
       'P3u', 'RDFV7', 'RDFC18', 'RDFV6', 'FNSA1', 'RDFC17', 'FNSA3', 'FNSA2', 'RDFC12', 'RDFC13', 'RDFC10', 'RDFC11',
       'P2p', 'RDFV4', 'MoRSEP22', 'RDFV3', 'MoRSEP18', 'RDFV2', 'RDFU21', 'RDFU20', 'RDFU23', 'RDFU22', 'RDFU25',
       'RDFM4', 'RDFU27', 'RDFU26', 'RDFU29', 'RDFU28', 'MoRSEN28', 'MoRSEN29', 'RDFV19', 'RDFV18', 'WPSA2', 'RDFV16',
       'RDFV15', 'RDFV14', 'RDFV13', 'RDFV12', 'RDFV11', 'RDFV10', 'MoRSEN26', 'MoRSEP23', 'MoRSEN27', 'MoRSEN24',
       'MoRSEP25', 'MoRSEN25', 'MoRSEP26', 'MoRSEE8', 'MoRSEE9', 'MoRSEE6', 'MoRSEN22', 'MoRSEE4', 'MoRSEP27',
       'MoRSEE2', 'MoRSEE3', 'RDFU24', 'MoRSEN20', 'MoRSEC9', 'ASPAN', 'RDFE10', 'MoRSEN21', 'Te', 'Vm', 'Vp',
       'MoRSEV18', 'PPSA3', 'Vv', 'RDFE13', 'E2u', 'RDFC30', 'E2v', 'P1m', 'MoRSEV12', 'MoRSEP15', 'MoRSEP17',
       'MoRSEU8', 'MoRSEU9', 'MoRSEU6', 'MoRSEU7', 'MoRSEU4', 'MoRSEU5', 'MoRSEU2', 'MoRSEU3', 'MoRSEU1', 'RDFM29',
       'RDFM28', 'MoRSEN6', 'RDFM21', 'RDFM20', 'RDFM23', 'RDFM22', 'RDFM25', 'RDFM24', 'RDFM27', 'RDFM26', 'RDFM2',
       'RDFM3', 'RDFV5', 'RDFM1', 'RDFM6', 'RDFM7', 'RDFV1', 'RDFM5', 'MoRSEP20', 'MoRSEP21', 'RDFM8', 'RDFM9',
       'MoRSEP24', 'RDFE8', 'RDFV9', 'RDFV8', 'MoRSEC8', 'RDFV17', 'RDFM17', 'WPSA3', 'AGDD', 'MoRSEC1', 'MoRSEC2',
       'MoRSEC3', 'MoRSEC4', 'MoRSEC5', 'MoRSEC6', 'MoRSEC7', 'MoRSEE29', 'MoRSEE28', 'WPSA1', 'MoRSEC29', 'MoRSEE21',
       'MoRSEE20', 'MoRSEE23', 'RDFM13', 'MoRSEE25', 'MoRSEE24', 'MoRSEE27', 'MoRSEE26', 'MoRSEC27', 'MoRSEC26', 'Ae',
       'RDFE1', 'RDFE6', 'RDFE7', 'RDFE4', 'MoRSEV28', 'MoRSEV29', 'MoRSEV26', 'MoRSEV27', 'MoRSEV24', 'MoRSEC20',
       'MoRSEV22', 'MoRSEV23', 'MoRSEV20', 'MoRSEV21', 'MoRSEC12', 'MoRSEC13', 'MoRSEC10', 'MoRSEC11', 'MoRSEC16',
       'MoRSEC17', 'MoRSEC14', 'MoRSEC15', 'P2m', 'MoRSEC18', 'MoRSEC19', 'RDFE30', 'RDFE21', 'RDFE20', 'RDFE23',
       'RDFE22', 'RDFE25', 'RDFE24', 'RDFE27', 'RDFE26', 'RDFE29', 'RDFE28', 'RDFP20', 'RDFP21', 'RDFP22', 'RDFP23',
       'RDFP24', 'RDFP25', 'RDFP26', 'RDFP27', 'RDFP28', 'RDFP29', 'MoRSEV6', 'MoRSEN13', 'MoRSEV30', 'Dv', 'RDFV26',
       'RDFV27', 'RDFV24', 'RDFV25', 'RDFV22', 'RDFV23', 'RDFV20', 'RDFV21', 'L2e', 'MoRSEV5', 'RDFV28', 'RDFV29',
       'MoRSEP3', 'P1u', 'rygr', 'Ve', 'MoRSEE7', 'MoRSEV4', 'MoRSEP2', 'FrTATP', 'MoRSEE5', 'P2v', 'ASA', 'MoRSEC21',
       'MoRSEV3', 'MoRSEE1', 'E3v', 'E3u', 'Ke', 'E3p', 'Km', 'MoRSEV7', 'E3e', 'Kp', 'Kv', 'Ku', 'MoRSEV2', 'RDFC8',
       'E3m', 'RDFC9', 'MoRSEM21', 'MoRSEM20', 'MoRSEM23', 'MoRSEM22', 'MoRSEM25', 'MoRSEM24', 'MoRSEM27', 'MoRSEM26',
       'MoRSEM29', 'MoRSEM28', 'MoRSEN19', 'MoRSEN18', 'L2m', 'MoRSEV9', 'MoRSEV8', 'MoRSEU29', 'MoRSEU28', 'L2u',
       'MoRSEV1', 'MoRSEN12', 'RDFC5', 'MoRSEU21', 'MoRSEU20', 'MoRSEU23', 'MoRSEU22', 'MoRSEU25', 'MoRSEU24',
       'MoRSEU27', 'MoRSEU26', 'RDFC7', 'MoRSEV25', 'Tv', 'Am', 'SEig', 'Tu', 'Tp', 'RDFC1', 'Tm', 'RDFC2', 'PNSA3',
       'PNSA2', 'PNSA1', 'RDFC3', 'MSA', 'MoRSEE22', 'L1m', 'De', 'MoRSEE10', 'MoRSEE11', 'MoRSEE12', 'MoRSEE13',
       'MoRSEP9', 'MoRSEE15', 'MoRSEE16', 'MoRSEE17', 'MoRSEE18', 'MoRSEE19', 'Du', 'Dp', 'Vu', 'P2e', 'E1p', 'E1u',
       'E1v', 'E1m', 'RNCS', 'MoRSEP28', 'MoRSEN30', 'E1e', 'MoRSEN7', 'W3D', 'RDFU18', 'RDFU19', 'RDFU14', 'RDFU15',
       'RDFU16', 'RDFU17', 'RDFU10', 'RDFU11', 'RDFU12', 'RDFU13', 'Ap', 'Au', 'RDFC4', 'MoRSEP29', 'Av', 'L3m',
       'RDFM30', 'MoRSEE14', 'E2m']



def getLdesc (typeDesc):

    lout = []
    if typeDesc == "1D2D":
        lout = list(constitution._constitutional.keys()) + list(molproperty._molProperty.keys()) + list(topology._topology.keys()) +\
                list(connectivity._connectivity.keys()) + list(kappa._kappa.keys()) + list(bcut._bcut.keys()) + list(basak._basak.keys()) +\
                list(EStateGlobal._EState.keys()) + list(moreaubroto._MBA.keys()) + list(moran._moran.keys()) + list(geary._geary.keys()) +\
                list(charge._charge.keys()) + list(moe._moe.keys()) + list(morgan._morgan.keys())

    elif typeDesc == "3D":
        lout = L3D

    return lout



class Descriptor:

    # mol have to be clean before
    def __init__(self, SMICLEAN, pdesc):
        self.smi = SMICLEAN
        self.mol = Chem.MolFromSmiles(SMICLEAN)
        self.err = 0
        self.pdesc = pdesc


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


    def computeAll2D(self):
        if path.exists(self.pdesc + "_2D.txt"):
            if path.getsize(self.pdesc + "_2D.txt") > 100:
                ddesc = toolbox.loadMatrixToDict(self.pdesc + "_2D.txt")
                self.all2D = ddesc
                return
            else:
                remove(self.pdesc + "_2D.txt")

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



    def computeAll3D(self,lcoordinates):

        if path.exists(self.pdesc + "_3D.txt"):
            if path.getsize(self.pdesc + "_3D.txt") > 100:
                self.all3D = toolbox.loadMatrixToDict(self.pdesc + "_3D.txt")
                return
            else:
                remove(self.pdesc + "_3D.txt")


        self.geo3D = GetGeo3D(lcoordinates)
        self.CPSA3D = GetCPSA3D(lcoordinates)
        self.rdf3D = GetRdf3D(lcoordinates)
        self.morse3D = GetMorse3D(lcoordinates)
        self.whim3D = GetWhim3D(lcoordinates)

        # combine 3D
        self.all3D = {}
        self.all3D.update(deepcopy(self.geo3D))
        self.all3D.update(deepcopy(self.CPSA3D))
        self.all3D.update(deepcopy(self.rdf3D))
        self.all3D.update(deepcopy(self.morse3D))
        self.all3D.update(deepcopy(self.whim3D))


    def writeMatrix(self, typedesc):
        if typedesc == "2D":
            if "all2D" in self.__dict__:
                filin = open(self.pdesc + "_2D.txt", "w")
                filin.write("%s\n"%("\t".join(self.all2D.keys())))
                filin.write("%s\n"%("\t".join([str(self.all2D[k]) for k in self.all2D.keys()])))
                filin.close()

        if typedesc == "3D":
            if "all3D" in self.__dict__:
                filin = open(self.pdesc + "_3D.txt", "w")
                filin.write("%s\n"%("\t".join(self.all3D.keys())))
                filin.write("%s\n"%("\t".join([str(self.all3D[k]) for k in self.all3D.keys()])))
                filin.close()



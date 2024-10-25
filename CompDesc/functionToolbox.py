from os import listdir, remove, makedirs, path, system, name
from shutil import rmtree
from re import search
from random import randint


#########Folder management#######
#################################
def cleanFolder(prin):

    lfiles = listdir(prin)
    if len(lfiles) != 0:
        for filin in lfiles:
            # problem with folder
            try: remove(prin + filin)
            except: rmtree(prin + filin)

    return prin


def createFolder(prin, clean=0):

    if not path.exists(prin):
        makedirs(prin)

    if clean == 1:
        cleanFolder(prin)

    return prin


###### LOADING MATRIX #######
#############################
def loadMatrixToDict(pmatrixIn, sep ="\t"):

    filin = open(pmatrixIn, "r")
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1)-1):
        lheaders.append("val")

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        #kin = lvalues[0]
        #dout[kin] = {}
        j = 0
        if len(lvalues) != len(lheaders):
            print(lineMat)
            print(llinesMat[i])
            print(lvalues)
            print("Error => nb element", i)
            print(len(lvalues))
            print(len(lheaders))

        jmax = len(lheaders)
        while j < jmax:
            dout[lheaders[j]] = lvalues[j]
            j += 1
        i += 1

    return dout


def formatLine(linein):

    linein = linein.replace("\n", "")
    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == '"' and flagchar == 0:
            flagchar = 1
        elif linein[i] == '"' and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + " "
        else:
            linenew = linenew + linein[i]
        i += 1

    linenew = linenew.replace('\"', "")
    return linenew

def parseSDFfor3DdescComputation(mol_block):
    """
    Read the coordinates and charge of each atom in molecule from .sdf file.
    """

    dchargeSDF = {7: -3.0, 6: -2.0, 5: -1.0, 0: 0.0, 3: 1.0, 2: 2.0, 1: 3.0}  # and 4 for radical

    latoms = []

    llines = str(mol_block).split("\n")

    # start at line 5 classical format
    for AtBlock in llines[4:]:
        if len(AtBlock) != 69 and len(AtBlock) != 51:
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

    return latoms

########RUN EXTERNAL SOFTWARE#########
######################################
def babelConvertMoltoSDF(pmolin, psdfout, update=1):

    if path.exists(psdfout) and update == 1:
        remove(psdfout)

    if not path.exists(psdfout):
        if name == "nt":
            print("TO DO FIX command line for window with openbabel 3.0.0 - l138 functionToolbox.py")
            cmd_convert = '"C:/Program Files (x86)/OpenBabel-2.3.1/babel.exe" ' + pmolin + " " + psdfout
            print(cmd_convert)
        else:
            cmd_convert = "obabel -imol " + pmolin + " -osdf -O " + psdfout + " 2>/dev/null"
            print(cmd_convert)
        system(cmd_convert)


####### PADEL and OPERA RUN #########
#####################################
PPADEL = "/usr/local/bin/OPERA/application/padel-full-1.00.jar"
PCDK = "/usr/local/bin/OPERA/application/CDKDescUI-2.0.jar"
OPERA = "/usr/local/bin/OPERA/application/run_OPERA.sh"
MATLAB = "/usr/local/MATLAB/MATLAB_Runtime/v99"

def runCDKDesc(psmi, prout, psoft):
    if psoft == "":
        psoft = PCDK
    a = randint(0, 100000000)
    pfilout = prout + str(a) + ".csv"
    if path.exists(pfilout) and path.getsize(pfilout) > 50:
        return pfilout

    if prout == "":
        return "ERROR - CDK2 Input"
    else:
        cmd = "java -jar %s -b -t all -o %s %s"%(psoft, pfilout, psmi)
        print(cmd)
        system(cmd)

    return pfilout


def runPadelDesc(prin, psoft):
    if psoft == "":
        psoft = PPADEL
    a = randint(0, 100000000)
    pfilout = prin + str(a) + ".csv"
    if path.exists(pfilout) and path.getsize(pfilout) > 50:
        return pfilout

    if prin == "":
        return "ERROR - Padel Input"
    else:
        cmd = "java -Djava.awt.headless=true -jar " + psoft + " -2d -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 10000 -dir " + str(
            prin) + " -file " + pfilout
        print(cmd)
        system(cmd)

    return pfilout

def runPadelFP(prin, psoft):#, pxml):
    if psoft == "":
        psoft = PPADEL

    # define pxml
    pxml = prin + "desc_fp.xml"
    fxml = open(pxml, "w")
    fxml.write("<Root>\n        <Group name=\"Fingerprint\">\n        <Descriptor name=\"Fingerprinter\" value=\"true\"/>\n        <Descriptor name=\"ExtendedFingerprinter\" value=\"true\"/>\n        <Descriptor name=\"EStateFingerprinter\" value=\"true\"/>\n        <Descriptor name=\"GraphOnlyFingerprinter\" value=\"true\"/>\n        <Descriptor name=\"MACCSFingerprinter\" value=\"true\"/>\n        <Descriptor name=\"PubchemFingerprinter\" value=\"true\"/>\n        <Descriptor name=\"SubstructureFingerprinter\" value=\"true\"/>\n        <Descriptor name=\"SubstructureFingerprintCount\" value=\"false\"/>\n        <Descriptor name=\"KlekotaRothFingerprinter\" value=\"true\"/>\n        <Descriptor name=\"KlekotaRothFingerprintCount\" value=\"false\"/>\n        <Descriptor name=\"AtomPairs2DFingerprinter\" value=\"true\"/>\n        <Descriptor name=\"AtomPairs2DFingerprintCount\" value=\"false\"/>\n    </Group>\n</Root>")
    fxml.close()
    # define a random name 
    a = randint(0, 100000000)
    pfilout = prin + str(a) + ".csv"
    if path.exists(pfilout) and path.getsize(pfilout) > 50:
        return pfilout

    if prin == "":
        return "ERROR - Padel Input"
    else:
        #pxml = "./doc/desc_fp.xml"
        #pxml = path.abspath(pxml)
        cmd = "java -Djava.awt.headless=true -jar %s -fingerprints -descriptortypes %s -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 10000 -dir %s -file %s" % (
        psoft, pxml, str(prin), pfilout)
        #cmd = "java -jar %s -fingerprints -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 10000 -dir %s -file %s" % (psoft, str(prin), pfilout)
        print(cmd)
        system(cmd)

    return pfilout


def runOPERA(p2Ddesc, pfp, pCDKdesc, pfilout, popera, pmatlab, onlyPhysChem=0, update = 0):
    """
    Compute only physico chemical properties
    """

    if popera == "":
        print("Use the default")
        popera = OPERA
    if pmatlab == "":
        pmatlab = MATLAB


    if path.exists(pfilout) and update == 0:
        return pfilout

    if onlyPhysChem == 1:
        cmd = "%s %s -d %s -fp %s -cdk %s -o %s -StrP -BCF -BP -logP -MP -VP -WS -AOH -BioDeg -RB -HL -KM -KOA -Koc -RT -logD -FuB -Clint -pKa"%(popera, pmatlab, p2Ddesc, pfp, pCDKdesc, pfilout)
    else:
        cmd = "%s %s -d %s -fp %s -cdk %s -o %s -StrP -BCF -BP -logP -MP -VP -WS -AOH -BioDeg -RB -HL -KM -KOA -Koc -RT -logD -FuB -Clint -pKa -CERAPP -CoMPARA -CATMoS"%(popera, pmatlab, p2Ddesc, pfp, pCDKdesc, pfilout)
  
    print (cmd)
    system(cmd)

    return pfilout

def runOPERAFromSmi(pSMI, pfilout, popera="", pmatlab="", update = 0):

    if popera == "":
        popera = OPERA
    if pmatlab == "":
        pmatlab = MATLAB


    if path.exists(pfilout) and update == 0:
        return pfilout

    cmd = "%s %s -s %s -o %s -a -st -v 1"%(popera, pmatlab, pSMI, pfilout)
    print (cmd)
    system(cmd)

    return pfilout   


def runOPERAFromChem(psmi, pfilout, popera, pmatlab, update=0):
    
    if popera == "":
        popera = OPERA
    if pmatlab == "":
        pmatlab = MATLAB

    if path.exists(pfilout) and update == 0:
        return pfilout

    cmd = "%s %s -s %s -o %s -a -StrP -BCF -BP -logP -MP -VP -WS -AOH -BioDeg -RB -HL -KM -KOA -Koc -RT -logD -FuB -Clint -pKa"%(popera, pmatlab, psmi, pfilout)


    system(cmd)

    return pfilout

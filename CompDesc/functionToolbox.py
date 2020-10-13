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

def parseSDFfor3DdescComputation(p3Dsdf):
    """
    Read the coordinates and charge of each atom in molecule from .sdf file.
    """

    dchargeSDF = {7: -3.0, 6: -2.0, 5: -1.0, 0: 0.0, 3: 1.0, 2: 2.0, 1: 3.0}  # and 4 for radical

    latoms = []

    filin = open(p3Dsdf, 'r')
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
OPERA = "/usr/local/bin/OPERA/application/run_OPERA.sh"
MATLAB = "/usr/local/MATLAB/MATLAB_Runtime/v94"
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
        cmd = "java -jar " + psoft + " -2d -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 10000 -dir " + str(
            prin) + " -file " + pfilout
        print(cmd)
        system(cmd)

    return pfilout

def runPadelFP(prin, psoft, pxml):
    if psoft == "":
        psoft = PPADEL

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
        cmd = "java -jar %s -fingerprints -descriptortypes %s -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 10000 -dir %s -file %s" % (
        psoft, pxml, str(prin), pfilout)
        print(cmd)
        system(cmd)

    return pfilout


def runOPERA(p2Ddesc, pfp, pfilout, popera, pmatlab, update = 0):

    if popera == "":
        popera = OPERA
    if pmatlab == "":
        pmatlab = MATLAB


    if path.exists(pfilout) and update == 0:
        return pfilout

    cmd = "%s %s -d %s -fp %s -o %s -StrP -BCF -BP -logP -MP -VP -WS -AOH -BioDeg -RB -HL -KM -KOA -Koc -RT -logD"%(popera, pmatlab, p2Ddesc, pfp, pfilout)
    print (cmd)
    system(cmd)

    return pfilout


from os import listdir, remove, makedirs, path, system
from shutil import rmtree
from re import search


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
def babelConvertMoltoSDF(pmolin, psdfout, window=0, update=1):

    if not path.exists(psdfout) or update==1:
        if window ==1:
            cmd_convert = '"C:/Program Files (x86)/OpenBabel-2.3.1/babel.exe" ' + pmolin + " " + psdfout
            print(cmd_convert)
        else:
            cmd_convert = "babel " + pmolin + " " + psdfout + " 2>/dev/null"
            print(cmd_convert)
        system(cmd_convert)

from rdkit import Chem
from .AtomProperty import GetRelativeAtomicProperty

import numpy


################################################################

def getMoranAutocorrelation(mol,lag=1,propertylabel='m'):
    """
    #################################################################
    Calculation of Moran autocorrelation descriptors based on
    different property weights.
    #################################################################
    """

    Natom = mol.GetNumAtoms()
    
    prolist = []
    for i in mol.GetAtoms():
        temp = GetRelativeAtomicProperty(i.GetSymbol(), propertyname=propertylabel)
        prolist.append(temp)
        
    aveweight=sum(prolist)/Natom
    
    tempp=[numpy.square(x-aveweight) for x in prolist]   
    
    GetDistanceMatrix=Chem.GetDistanceMatrix(mol)
    res = 0.0
    index = 0
    for i in range(Natom):
        for j in range(Natom):  
            if GetDistanceMatrix[i,j] == lag:
                atom1 = mol.GetAtomWithIdx(i)
                atom2 = mol.GetAtomWithIdx(j)
                temp1 = GetRelativeAtomicProperty(element=atom1.GetSymbol(), propertyname=propertylabel)
                temp2 = GetRelativeAtomicProperty(element=atom2.GetSymbol(), propertyname=propertylabel)
                res = res + (temp1-aveweight)*(temp2-aveweight)
                index = index+1
            else:
                res = res+0.0
                
                
    if sum(tempp) == 0 or index == 0:
        result = 0.0
    else:
        result = (res/index)/(sum(tempp)/Natom)
                
    return round(result, 6)

#############################

def getMATSm(mol):
    dres = {}
    for i in range(8):
        dres['MATSm'+str(i+1)] = getMoranAutocorrelation(mol, lag=i+1, propertylabel='m')
    return dres

def getMATSv(mol):
    dres = {}
    for i in range(8):
        dres['MATSv'+str(i+1)] = getMoranAutocorrelation(mol, lag=i+1, propertylabel='V')
    return dres

def getMATSe(mol):
    dres = {}
    for i in range(8):
        dres['MATSe'+str(i+1)] = getMoranAutocorrelation(mol, lag=i+1, propertylabel='En')
    return dres

def getMATSp(mol):
    dres = {}
    for i in range(8):
        dres['MATSp'+str(i+1)] = getMoranAutocorrelation(mol, lag=i+1, propertylabel='alapha')
    return dres


_moran = {"MATSm1":getMATSm,
          "MATSm2":getMATSm,
          "MATSm3":getMATSm,
          "MATSm4":getMATSm,
          "MATSm5":getMATSm,
          "MATSm6":getMATSm,
          "MATSm7":getMATSm,
          "MATSm8":getMATSm,
          "MATSv1":getMATSv,
          "MATSv2":getMATSv,
          "MATSv3":getMATSv,
          "MATSv4":getMATSv,
          "MATSv5":getMATSv,
          "MATSv6":getMATSv,
          "MATSv7":getMATSv,
          "MATSv8":getMATSv,
          "MATSe1":getMATSe,
          "MATSe2":getMATSe,
          "MATSe3":getMATSe,
          "MATSe4":getMATSe,
          "MATSe5":getMATSe,
          "MATSe6":getMATSe,
          "MATSe7":getMATSe,
          "MATSe8":getMATSe,
          "MATSp1":getMATSp,
          "MATSp2":getMATSp,
          "MATSp3":getMATSp,
          "MATSp4":getMATSp,
          "MATSp5":getMATSp,
          "MATSp6":getMATSp,
          "MATSp7":getMATSp,
          "MATSp8":getMATSp}


def GetMATS(mol):
    dres={}
    dres.update(getMATSm(mol))
    dres.update(getMATSv(mol))
    dres.update(getMATSe(mol))
    dres.update(getMATSp(mol))
    return dres

from rdkit import Chem
from .AtomProperty import GetRelativeAtomicProperty
import numpy

def getGearyAutocorrelation(mol, lag=1, propertylabel='m'):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on
    different property weights.
    #################################################################
    """

    Natom = mol.GetNumAtoms()
    prolist = []
    for i in mol.GetAtoms():
        temp = GetRelativeAtomicProperty(i.GetSymbol(), propertyname=propertylabel)
        prolist.append(temp)
        
    aveweight=sum(prolist)/Natom
    
    tempp = [numpy.square(x-aveweight) for x in prolist]
    
    GetDistanceMatrix = Chem.GetDistanceMatrix(mol)
    res = 0.0
    index = 0
    for i in range(Natom):
        for j in range(Natom):  
            if GetDistanceMatrix[i,j] == lag:
                atom1 = mol.GetAtomWithIdx(i)
                atom2 = mol.GetAtomWithIdx(j)
                temp1 = GetRelativeAtomicProperty(element=atom1.GetSymbol(), propertyname=propertylabel)
                temp2 = GetRelativeAtomicProperty(element=atom2.GetSymbol(), propertyname=propertylabel)
                res = res+numpy.square(temp1-temp2)
                index = index+1
            else:
                res = res + 0.0
                
                
    if sum(tempp)==0 or index==0:
        result=0
    else:
        result=(res/index/2)/(sum(tempp)/(Natom-1))
                
    return round(result,6)


###########

def getGATSm(mol):
    dres={}
    for i in range(8):
        dres['GATSm'+str(i+1)] = getGearyAutocorrelation(mol, lag=i+1, propertylabel='m')
    return dres

def getGATSv(mol):
    dres={}
    for i in range(8):
        dres['GATSv'+str(i+1)] = getGearyAutocorrelation(mol, lag=i+1, propertylabel='V')
    return dres

def getGATSe(mol):
    dres={}
    for i in range(8):
        dres['GATSe'+str(i+1)] = getGearyAutocorrelation(mol, lag=i+1, propertylabel='En')
    return dres

def getGATSp(mol):
    dres={}
    for i in range(8):
        dres['GATSp'+str(i+1)] = getGearyAutocorrelation(mol, lag=i+1, propertylabel='alapha')
    return dres


_geary = {"GATSm1":getGATSm,
          "GATSm2":getGATSm,
          "GATSm3":getGATSm,
          "GATSm4":getGATSm,
          "GATSm5":getGATSm,
          "GATSm6":getGATSm,
          "GATSm7":getGATSm,
          "GATSm8":getGATSm,
          "GATSv1":getGATSv,
          "GATSv2":getGATSv,
          "GATSv3":getGATSv,
          "GATSv4":getGATSv,
          "GATSv5":getGATSv,
          "GATSv6":getGATSv,
          "GATSv7":getGATSv,
          "GATSv8":getGATSv,
          "GATSe1":getGATSe,
          "GATSe2":getGATSe,
          "GATSe3":getGATSe,
          "GATSe4":getGATSe,
          "GATSe5":getGATSe,
          "GATSe6":getGATSe,
          "GATSe7":getGATSe,
          "GATSe8":getGATSe,
          "GATSp1":getGATSp,
          "GATSp2":getGATSp,
          "GATSp3":getGATSp,
          "GATSp4":getGATSp,
          "GATSp5":getGATSp,
          "GATSp6":getGATSp,
          "GATSp7":getGATSp,
          "GATSp8":getGATSp}


def GetGATS(mol):
    res={}
    res.update(getGATSm(mol))
    res.update(getGATSv(mol))
    res.update(getGATSe(mol))
    res.update(getGATSp(mol))
    return res

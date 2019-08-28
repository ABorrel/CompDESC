from rdkit import Chem
from .AtomProperty import GetRelativeAtomicProperty

import numpy


def getMBautocorelation(mol, lag=1, propertylabel='m'):
    """
    #################################################################
    Calculation of Moreau-Broto autocorrelation descriptors based on 
    different properties (mass, volume, electro, polarizability).
    #################################################################
    """

    Natom = mol.GetNumAtoms()
    GetDistanceMatrix = Chem.GetDistanceMatrix(mol)
    res=0.0
    for i in range(Natom):
        for j in range(Natom):  
            if GetDistanceMatrix[i,j] == lag:
                atom1=mol.GetAtomWithIdx(i)
                atom2=mol.GetAtomWithIdx(j)
                temp1=GetRelativeAtomicProperty(element=atom1.GetSymbol(),propertyname=propertylabel)
                temp2=GetRelativeAtomicProperty(element=atom2.GetSymbol(),propertyname=propertylabel)
                res=res+temp1*temp2
            else:
                res=res+0.0
                
    return round(numpy.log(res/2+1),6)

###############

def getMBAm(mol):
    """
    #################################################################
    Calculation of Moreau-Broto autocorrelation descriptors based on 
    carbon-scaled atomic mass.
    #################################################################
    """
    res={}
    for i in range(8):
        res['MBAm'+str(i+1)]=getMBautocorelation(mol,lag=i+1,propertylabel='m')
    return res


def getMBAv(mol):
    """
    #################################################################
    Calculation of Moreau-Broto autocorrelation descriptors based on 
    carbon-scaled atomic van der Waals volume.
    #################################################################
    """
    res={}
    for i in range(8):
        res['MBAv'+str(i+1)]=getMBautocorelation(mol,lag=i+1,propertylabel='V')
    return res

def getMBAe(mol):
    """
    #################################################################
    Calculation of Moreau-Broto autocorrelation descriptors based on 
    carbon-scaled atomic Sanderson electronegativity.
    #################################################################
    """
    res={}
    for i in range(8):
        res['MBAe'+str(i+1)]=getMBautocorelation(mol, lag=i+1, propertylabel='En')
    return res

def getMBAp(mol):
    """
    #################################################################
    Calculation of Moreau-Broto autocorrelation descriptors based on 
    carbon-scaled atomic polarizability.
    #################################################################
    """
    res={}
    for i in range(8):
        res['MBAp'+str(i+1)]=getMBautocorelation(mol,lag=i+1,propertylabel='alapha')
    return res


_MBA = {"MBAm1":getMBAm,
        "MBAm2":getMBAm,
        "MBAm3":getMBAm,
        "MBAm4":getMBAm,
        "MBAm5":getMBAm,
        "MBAm6":getMBAm,
        "MBAm7":getMBAm,
        "MBAm8":getMBAm,
        "MBAv1":getMBAv,
        "MBAv2":getMBAv,
        "MBAv3":getMBAv,
        "MBAv4":getMBAv,
        "MBAv5":getMBAv,
        "MBAv6":getMBAv,
        "MBAv7":getMBAv,
        "MBAv8":getMBAv,
        "MBAe1":getMBAe,
        "MBAe2":getMBAe,
        "MBAe3":getMBAe,
        "MBAe4":getMBAe,
        "MBAe5":getMBAe,
        "MBAe6":getMBAe,
        "MBAe7":getMBAe,
        "MBAe8":getMBAe,
        "MBAp1":getMBAp,
        "MBAp2":getMBAp,
        "MBAp3":getMBAp,
        "MBAp4":getMBAp,
        "MBAp5":getMBAp,
        "MBAp6":getMBAp,
        "MBAp7":getMBAp,
        "MBAp8":getMBAp}


def GetMBA(mol):
    dres={}
    dres.update(getMBAm(mol))
    dres.update(getMBAv(mol))
    dres.update(getMBAe(mol))
    dres.update(getMBAp(mol))
    return dres

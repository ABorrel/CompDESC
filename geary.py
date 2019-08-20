from rdkit import Chem
from .AtomProperty import GetRelativeAtomicProperty

import numpy

Version=1.0


################################################################
def _CalculateGearyAutocorrelation(mol,lag=1,propertylabel='m'):
    """
    #################################################################
    **Internal used only**
    
    Calculation of Geary autocorrelation descriptors based on 
    
    different property weights.
    
    Usage:
        
    res=_CalculateGearyAutocorrelation(mol,lag=1,propertylabel='m')
    
    Input: mol is a molecule object.
    
    lag is the topological distance between atom i and atom j.
    
    propertylabel is the weighted property.
    
    Output: res is a numeric value.
    #################################################################
    """

    Natom=mol.GetNumAtoms()
    
    prolist=[]
    for i in mol.GetAtoms():
        temp=GetRelativeAtomicProperty(i.GetSymbol(),propertyname=propertylabel)
        prolist.append(temp)
        
    aveweight=sum(prolist)/Natom
    
    tempp=[numpy.square(x-aveweight) for x in prolist]   
    
    GetDistanceMatrix=Chem.GetDistanceMatrix(mol)
    res=0.0
    index=0
    for i in range(Natom):
        for j in range(Natom):  
            if GetDistanceMatrix[i,j]==lag:
                atom1=mol.GetAtomWithIdx(i)
                atom2=mol.GetAtomWithIdx(j)
                temp1=GetRelativeAtomicProperty(element=atom1.GetSymbol(),propertyname=propertylabel)
                temp2=GetRelativeAtomicProperty(element=atom2.GetSymbol(),propertyname=propertylabel)
                res=res+numpy.square(temp1-temp2)
                index=index+1
            else:
                res=res+0.0
                
                
    if sum(tempp)==0 or index==0:
        result=0
    else:
        result=(res/index/2)/(sum(tempp)/(Natom-1))
                
    return round(result,3)


def CalculateGearyAutoMass(mol):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on 
    
    carbon-scaled atomic mass.
    
    Usage:
    
    res=CalculateMoranAutoMass(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['GATSm'+str(i+1)]=_CalculateGearyAutocorrelation(mol,lag=i+1,propertylabel='m')
    
    
    return res


def CalculateGearyAutoVolume(mol):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on 
    
    carbon-scaled atomic van der Waals volume.

    Usage:
    
    res=CalculateGearyAutoVolume(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['GATSv'+str(i+1)]=_CalculateGearyAutocorrelation(mol,lag=i+1,propertylabel='V')
    
    
    return res

def CalculateGearyAutoElectronegativity(mol):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on 
    
    carbon-scaled atomic Sanderson electronegativity.
    
    Usage:
    
    res=CalculateGearyAutoElectronegativity(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['GATSe'+str(i+1)]=_CalculateGearyAutocorrelation(mol,lag=i+1,propertylabel='En')
    
    
    return res

def CalculateGearyAutoPolarizability(mol):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on 
    
    carbon-scaled atomic polarizability.
    
    Usage:
    
    res=CalculateGearyAutoPolarizability(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['GATSp'+str(i+1)]=_CalculateGearyAutocorrelation(mol,lag=i+1,propertylabel='alapha')
    
    
    return res


def GetGearyAuto(mol):
    """
    #################################################################
    Calcualate all Geary autocorrelation descriptors.

    (carbon-scaled atomic mass, carbon-scaled atomic van der Waals volume,
     
    carbon-scaled atomic Sanderson electronegativity,
     
    carbon-scaled atomic polarizability)
    
    Usage:
    
    res=GetGearyAuto(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing all geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    res.update(CalculateGearyAutoMass(mol))
    res.update(CalculateGearyAutoVolume(mol))
    res.update(CalculateGearyAutoElectronegativity(mol))
    res.update(CalculateGearyAutoPolarizability(mol))
    
    return res

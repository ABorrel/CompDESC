from rdkit import Chem
from .AtomProperty import GetRelativeAtomicProperty
import numpy
import numpy.linalg


def getBurdenMatrix(mol, propertylabel='m'):
    """
    #################################################################
    Calculate Burden matrix and their eigenvalues.
    #################################################################
    """
    mol=Chem.AddHs(mol)
    Natom=mol.GetNumAtoms()

    AdMatrix=Chem.GetAdjacencyMatrix(mol)
    bondindex=numpy.argwhere(AdMatrix)
    AdMatrix1=numpy.array(AdMatrix,dtype=float)
    
    #The diagonal elements of B, Bii, are either given by 
    #the carbon normalized atomic mass,
    #van der Waals volume, Sanderson electronegativity,
    #and polarizability of atom i.

    for i in range(Natom):
        atom=mol.GetAtomWithIdx(i)
        temp=GetRelativeAtomicProperty(element=atom.GetSymbol(),propertyname=propertylabel)
        AdMatrix1[i,i]=round(temp,3)
        
    #The element of B connecting atoms i and j, Bij, 
    #is equal to the square root of the bond
    #order between atoms i and j.
    
    for i in bondindex:
        bond=mol.GetBondBetweenAtoms(int(i[0]),int(i[1]))
        if bond.GetBondType().name =='SINGLE':
            AdMatrix1[i[0],i[1]]=round(numpy.sqrt(1),3)
        if bond.GetBondType().name=="DOUBLE":
            AdMatrix1[i[0],i[1]]=round(numpy.sqrt(2),3)
        if bond.GetBondType().name=="TRIPLE":
            AdMatrix1[i[0],i[1]]=round(numpy.sqrt(3),3)
        if bond.GetBondType().name=="AROMATIC":
            AdMatrix1[i[0],i[1]]=round(numpy.sqrt(1.5),3)
    
    ##All other elements of B (corresponding non bonded 
    #atom pairs) are set to 0.001       
    bondnonindex=numpy.argwhere(AdMatrix==0)
    
    for i in bondnonindex:
        if i[0]!=i[1]:
            AdMatrix1[i[0],i[1]]=0.001
    return numpy.real(numpy.linalg.eigvals(AdMatrix1))
    

########################

def getbcutm(mol):
    """
    #################################################################
    Calculate Burden descriptors based on atomic mass.
    res--->dict type with 16 descriptors
    #################################################################
    """
    temp=getBurdenMatrix(mol,propertylabel='m')
    temp1=numpy.sort(temp[temp>=0])
    temp2=numpy.sort(numpy.abs(temp[temp<0]))
    
    if len(temp1)<8:
        temp1=numpy.concatenate((numpy.zeros(8),temp1))
    if len(temp2)<8:
        temp2=numpy.concatenate((numpy.zeros(8),temp2))
    
    bcut=["bcutm16","bcutm15","bcutm14","bcutm13","bcutm12","bcutm11","bcutm10",
          "bcutm9","bcutm8","bcutm7","bcutm6","bcutm5","bcutm4","bcutm3",
          "bcutm2","bcutm1"]
    bcutvalue=numpy.concatenate((temp2[-8:],temp1[-8:]))
    
    bcutvalue=[round(i,6) for i in bcutvalue]
    res=dict(zip(bcut,bcutvalue))
    return res

        

def getbcutv(mol):
    """
    #################################################################
    Calculate Burden descriptors based on atomic vloumes
    res-->dict type with 16 descriptors
    #################################################################
    """
    temp=getBurdenMatrix(mol,propertylabel='V')
    temp1=numpy.sort(temp[temp>=0])
    temp2=numpy.sort(numpy.abs(temp[temp<0]))
    
    if len(temp1)<8:
        temp1=numpy.concatenate((numpy.zeros(8),temp1))
    if len(temp2)<8:
        temp2=numpy.concatenate((numpy.zeros(8),temp2))
    
    bcut=["bcutv16","bcutv15","bcutv14","bcutv13","bcutv12","bcutv11","bcutv10",
          "bcutv9","bcutv8","bcutv7","bcutv6","bcutv5","bcutv4","bcutv3",
          "bcutv2","bcutv1"]
    bcutvalue=numpy.concatenate((temp2[-8:],temp1[-8:]))
    
    bcutvalue=[round(i,3) for i in bcutvalue]
    res=dict(zip(bcut,bcutvalue))
    return res



def getbcute(mol):
    """
    #################################################################
    Calculate Burden descriptors based on atomic electronegativity.
    res-->dict type with 16 descriptors
    #################################################################
    """
    temp=getBurdenMatrix(mol,propertylabel='En')
    temp1=numpy.sort(temp[temp>=0])
    temp2=numpy.sort(numpy.abs(temp[temp<0]))
    
    if len(temp1)<8:
        temp1=numpy.concatenate((numpy.zeros(8),temp1))
    if len(temp2)<8:
        temp2=numpy.concatenate((numpy.zeros(8),temp2))
    
    bcut=["bcute16","bcute15","bcute14","bcute13","bcute12","bcute11","bcute10",
          "bcute9","bcute8","bcute7","bcute6","bcute5","bcute4","bcute3",
          "bcute2","bcute1"]
    bcutvalue=numpy.concatenate((temp2[-8:],temp1[-8:]))
    
    bcutvalue=[round(i,3) for i in bcutvalue]
    res=dict(zip(bcut,bcutvalue))
    return res


def getbcutp(mol):
    """
    #################################################################
    Calculate Burden descriptors based on polarizability.
    res-->dict type with 16 descriptors
    #################################################################
    """
    temp=getBurdenMatrix(mol,propertylabel='alapha')
    temp1=numpy.sort(temp[temp>=0])
    temp2=numpy.sort(numpy.abs(temp[temp<0]))
    
    if len(temp1)<8:
        temp1=numpy.concatenate((numpy.zeros(8),temp1))
    if len(temp2)<8:
        temp2=numpy.concatenate((numpy.zeros(8),temp2))
    
    bcut=["bcutp16","bcutp15","bcutp14","bcutp13","bcutp12","bcutp11","bcutp10",
          "bcutp9","bcutp8","bcutp7","bcutp6","bcutp5","bcutp4","bcutp3",
          "bcutp2","bcutp1"]
    bcutvalue=numpy.concatenate((temp2[-8:],temp1[-8:]))
    
    bcutvalue=[round(i,3) for i in bcutvalue]
    res=dict(zip(bcut,bcutvalue))
    return res


_bcut = {"bcutm1": getbcutm,
         "bcutm2": getbcutm,
         "bcutm3": getbcutm,
         "bcutm4": getbcutm,
         "bcutm5": getbcutm,
         "bcutm6": getbcutm,
         "bcutm7": getbcutm,
         "bcutm8": getbcutm,
         "bcutm9": getbcutm,
         "bcutm10": getbcutm,
         "bcutm11": getbcutm,
         "bcutm12": getbcutm,
         "bcutm13": getbcutm,
         "bcutm14": getbcutm,
         "bcutm15": getbcutm,
         "bcutm16": getbcutm,
         "bcutv1": getbcutv,
         "bcutv2": getbcutv,
         "bcutv3": getbcutv,
         "bcutv4": getbcutv,
         "bcutv5": getbcutv,
         "bcutv6": getbcutv,
         "bcutv7": getbcutv,
         "bcutv8": getbcutv,
         "bcutv9": getbcutv,
         "bcutv10": getbcutv,
         "bcutv11": getbcutv,
         "bcutv12": getbcutv,
         "bcutv13": getbcutv,
         "bcutv14": getbcutv,
         "bcutv15": getbcutv,
         "bcutv16": getbcutv,
         "bcute1": getbcute,
         "bcute2": getbcute,
         "bcute3": getbcute,
         "bcute4": getbcute,
         "bcute5": getbcute,
         "bcute6": getbcute,
         "bcute7": getbcute,
         "bcute8": getbcute,
         "bcute9": getbcute,
         "bcute10": getbcute,
         "bcute11": getbcute,
         "bcute12": getbcute,
         "bcute13": getbcute,
         "bcute14": getbcute,
         "bcute15": getbcute,
         "bcute16": getbcute,
         "bcutp1": getbcutp,
         "bcutp2": getbcutp,
         "bcutp3": getbcutp,
         "bcutp4": getbcutp,
         "bcutp5": getbcutp,
         "bcutp6": getbcutp,
         "bcutp7": getbcutp,
         "bcutp8": getbcutp,
         "bcutp9": getbcutp,
         "bcutp10": getbcutp,
         "bcutp11": getbcutp,
         "bcutp12": getbcutp,
         "bcutp13": getbcutp,
         "bcutp14": getbcutp,
         "bcutp15": getbcutp,
         "bcutp16": getbcutp
         }


def GetBcut(mol):
    """
    #################################################################
    Calculate all 64 Burden descriptors
    res-->dict type
    #################################################################
    """
    dresult={}
    dresult.update(getbcutm(mol))
    dresult.update(getbcutv(mol))
    dresult.update(getbcute(mol))
    dresult.update(getbcutp(mol))
    return dresult




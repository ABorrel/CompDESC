from rdkit.Chem import rdchem
periodicTable = rdchem.GetPeriodicTable()
from math import pi

"""
Z: atomic number
L: principal quantum number -> row in the peiodic table
Zv: number of valence electrons
Rv: van der Waals atomic radius
Rc: covalent radius
m: atomic mass
V: van der Waals vloume
En: Sanderson electronegativity
alapha: atomic polarizability (10e-30 m3) -> from https://www.sciencedirect.com/bookseries/advances-in-atomic-and-molecular-physics/vol/13
IP: ionization potential (eV) -> ionization energy (http://www.knowledgedoor.com/2/elements_handbook/hydrogen.html)
EA: electron affinity (eV)
"""

################################################################################
AtomProperty={
    'H':{'L':1,'En':2.59,'alapha':0.67,'IP':13.598,'EA':0.754},
    'He':{'L':1,'En': 0.0,'alapha':0.20, 'IP': 24.59, 'EA':0.0},
    'Li':{'L':2,'En':0.89,'alapha':24.3,'IP':5.392,'EA':0.618},
    'Be':{'L':2,'En':1.81,'alapha':5.60,'IP':9.323,'EA':0.0},
    'B':{'L':2,'En':2.28,'alapha':3.03,'IP':8.298,'EA':0.277},
    'C':{'L':2,'En':2.75,'alapha':1.76,'IP':11.260,'EA':1.263},
    'N':{'L':2,'En':3.19,'alapha':1.10,'IP':14.534,'EA':0.0},
    'O':{'L':2,'En':3.65,'alapha':0.80,'IP':13.618,'EA':1.461},
    'F':{'L':2,'En':4.00,'alapha':0.56,'IP':17.423,'EA':3.401},
    'Ne':{'L':2,'En':4.50,'alapha':0.40, 'IP': 21.56, 'EA':0.0},
    'Na':{'L':3,'En':0.56,'alapha':23.6,'IP':5.139,'EA':0.548},
    'Mg':{'L':3,'En':1.32,'alapha':10.6,'IP':7.646,'EA':0.0},
    'Al':{'L':3,'En':1.71,'alapha':8.34,'IP':5.986,'EA':0.441},
    'Si':{'L':3,'En':2.14,'alapha':5.38,'IP':8.152,'EA':1.385},
    'P':{'L':3,'En':2.52,'alapha':3.63,'IP':10.487,'EA':0.747},
    'S':{'L':3,'En':2.96,'alapha':2.90,'IP':10.360,'EA':2.077},
    'Cl':{'L':3,'En':3.48,'alapha':2.18,'IP':12.968,'EA':3.613},
    'Ar':{'L':3,'En':3.31,'alapha':1.64,'IP':15.76, 'EA':0.0},
    'K':{'L':4,'En':0.45,'alapha':43.4,'IP':4.341,'EA':0.501},
    'Ca':{'L':4,'En':0.95,'alapha':25.0,'IP':6.113,'EA':0.018},
    'Sc':{'L':4,'En':1.02, 'alapha':16.9,'IP':6.56,'EA':0.188},
    'Ti':{'L':4,'En':1.50,'alapha':13.6, 'IP':6.83, 'EA':0.084},
    'V':{'L':4,'En':2.51, 'alapha':11.4,'IP': 6.75,'EA':0.525},
    'Cr':{'L':4,'En':1.66,'alapha':6.8,'IP':6.767,'EA':0.666},
    'Mn':{'L':4,'En':2.20,'alapha':8.6,'IP':7.434,'EA':0.0},
    'Fe':{'L':4,'En':2.20,'alapha':7.50,'IP':7.902,'EA':1.151},
    'Co':{'L':4,'En':2.56,'alapha':6.80,'IP':7.881,'EA':0.662},
    'Ni':{'L':4,'En':1.94,'alapha':6.5,'IP':7.640,'EA':1.156},
    'Cu':{'L':4,'En':1.95,'alapha':6.10,'IP':7.723,'EA':1.235},
    'Zn':{'L':4,'En':2.23,'alapha':7.08,'IP':9.394,'EA':0.0},
    'Ga':{'L':4,'En':2.42,'alapha':8.12,'IP':5.999,'EA':0.300},
    'Ge':{'L':4,'En':2.62,'alapha':6.07,'IP':7.900,'EA':1.233},
    'As':{'L':4,'En':2.82,'alapha':4.31,'IP':9.815,'EA':0.810},
    'Se':{'L':4,'En':3.01,'alapha':3.77,'IP':9.752,'EA':2.021},
    'Br':{'L':4,'En':3.22,'alapha':3.05,'IP':11.814,'EA':3.364},
    'Kr':{'L':4,'En':2.91,'alapha':2.48, 'IP':14.00,'EA':0.0},
    'Rb':{'L':5,'En':0.31,'alapha':47.3,'IP':4.177,'EA':0.486},
    'Sr':{'L':5,'En':0.72,'alapha':27.6,'IP':5.695,'EA':0.110},
    'Y':{'L':5,'En':0.65,'alapha':22.0,'IP':6.22,'EA':0.307},
    'Zr':{'L':5,'En':0.90,'alapha':18.0,'IP':6.63, 'EA':0.426},
    'Nb':{'L':5,'En':1.42,'alapha':14.0,'IP':6.76,'EA':0.893},
    'Mo':{'L':5,'En':1.15,'alapha':13.0,'IP':7.092,'EA':0.746},
    'Tc':{'L':5,'En':0.0,'alapha':10.0,'IP':7.28,'EA':0.55},
    'Ru':{'L':5,'En':0.0,'alapha':8.6,'IP':7.36,'EA':1.05},
    "Rh":{'L':5,'En':0.0,'alapha':7.6,'IP':7.46, 'EA':1.14},
    'Pd':{'L':5,'En':1.58,'alapha': 6.90,'IP':8.34,'EA':0.56},
    'Ag':{'L':5,'En':1.83,'alapha':6.30,'IP':7.576,'EA':1.302},
    'Cd':{'L':5,'En':1.98,'alapha':6.0,'IP':8.994,'EA':0.0},
    'In':{'L':5,'En':2.14,'alapha':4.5,'IP':5.786,'EA':0.300},
    'Sn':{'L':5,'En':2.30,'alapha':4.4,'IP':7.344,'EA':1.112},
    'Sb':{'L':5,'En':2.46,'alapha':4.0,'IP':8.64,'EA':1.07},
    'Te':{'L':5,'En':2.62,'alapha':3.9,'IP':9.01,'EA':1.971},
    'I':{'L':5,'En':2.78,'alapha':3.9,'IP':10.451,'EA':3.059},
    'Xe':{'L':5,'En':2.34,'alapha':4.04,'IP':12.13,'EA':0.0},
    'Cs':{'L':6,'En':0.22,'alapha':59.6, 'IP':3.89, 'EA':0.472},
    'Ba':{'L':6,'En':0.683,'alapha':39.7,'IP':5.21, 'EA':0.145},
    'La':{'L':6,'En':0.0,'alapha':37.0,'IP':5.58, 'EA':0.47},
    'Ce':{'L':6,'En':0.0,'alapha':36.0, 'IP':5.54, 'EA':0.955},
    #'Pr':{'L':6,'En':,'alapha':34.0},
    'Nd':{'L':6,'En':0.0,'alapha':32.0,'IP':5.525,'EA':0.518 },
    #'Pm':{'L':6,'En':,'alapha':30.0},
    'Sm':{'L':6,'En':0.0,'alapha':29.0, 'IP': 5.6437, 'EA': 0.518},
    'Eu':{'L':6,'En':0.0,'alapha':27.0, 'IP': 5.67, 'EA':0.864},
    'Gd':{'L':6,'En':2.00,'alapha':26.0,'IP':6.15,'EA':0.50},
    'Bk':{'L':6,'En':0.00,'alapha':39.0,'IP':6.198,'EA':0.0},
    #'Tb':{'alapha':25.0},
    'Dy':{'L':6, 'En':0.0,'alapha':25.0, 'IP': 5.94, 'EA': 0.518},
    #'Ho':{'alapha':23.0},
    #'Er':{'alapha':23.0},
    #'Tm':{22.0},
    #'Yb':{'alapha':22.0},
    #'Lu':{'alapha':20.0},
    #'Hf':{'alapha':15.0},
    #'Ta':{'alapha':13.0},
    #'W':{'alapha':10.0},
    #'Re':{'alapha':9.0},
    #'Os':{'alapha':8.0},
    #'Ir':{'alapha':7.0},
    'Pt':{'L':6,'Zv':10,'Rv':1.75,'Rc':1.28,'m':195.08,'V':22.45,'En':2.28,'alapha':6.30,'IP':9.00,'EA':2.128},
    'Au':{'L':6,'Zv':11,'Rv':1.66,'Rc':1.44,'m':196.97,'V':19.16,'En':2.65,'alapha':5.70,'IP':9.226,'EA':2.309},
    'Hg':{'L':6,'Zv':12,'Rv':1.55,'Rc':1.49,'m':200.59,'V':15.60,'En':2.20,'alapha':5.10,'IP':10.438,'EA':0.0},
    'Tl':{'L':6,'Zv':3,'Rv':1.96,'Rc':1.48,'m':204.38,'V':31.54,'En':2.25,'alapha':3.50,'IP':6.108,'EA':0.200},
    'Pb':{'L':6,'Zv':4,'Rv':2.02,'Rc':1.47,'m':207.20,'V':34.53,'En':2.29,'alapha':3.70,'IP':7.417,'EA':0.364},
    'Bi':{'L':6,'En':2.34,'alapha':4.0,'IP':7.289,'EA':0.946},
    'At':{'L':6,'En':0.0,'Zv':7,'Rv':2.88,'Rc':1.50,'m':85,'alapha':8.3,'IP':9.32,'EA':2.8}, # to double check for alpha from https://core.ac.uk/download/pdf/162160244.pdf
    'Ra':{'L':7,'En':0.0,'alapha':46.0,'IP':5.2784,'EA':0.10},
    'Th':{'L':7,'En':0.0,'alapha':50.0,'IP':6.3067,'EA':0.0},
}



def GetAbsoluteAtomicProperty(element='C',propertyname='m'):
    """
    Get the absolute property value with propertyname for the given atom.
    """

    if propertyname == "m":
        return periodicTable.GetAtomicWeight(element)
    elif propertyname == "V":
        r = periodicTable.GetRvdw(element)
        V = 4/3*pi*r**3
        return V
    elif propertyname == "Z":
        return periodicTable.GetAtomicNumber(element)
    elif propertyname == "Rv":
        return periodicTable.GetRvdw(element)
    elif propertyname == "Rc":
        return periodicTable.GetRb0(element)
    elif propertyname == "Zv":
        return periodicTable.GetDefaultValence(element)
    else:
        PropertyDic = AtomProperty[element]
        return PropertyDic[propertyname]



def GetRelativeAtomicProperty(element='C',propertyname='m'):
    """
    Get the absolute property value with propertyname for the given atom.
    """
    
    CpropertyDic = float(GetAbsoluteAtomicProperty('C', propertyname))
    PropertyDic = float(GetAbsoluteAtomicProperty(element, propertyname))
    
    return PropertyDic/CpropertyDic


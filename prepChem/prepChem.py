from re import search
from rdkit import Chem
import urllib.request
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

# frag not considered
# not use because only metal and add function is_metal
LSMILESREMOVE=["[C-]#N", "[Al+3]", "[Gd+3]", "[Pt+2]", "[Au+3]", "[Bi+3]", "[Al]", "[Si+4]", "[Fe]", "[Zn]", "[Fe+2]",
               "[Ru+8]", "[Fe+]", "[Sr++]", "[Fe+3]", "[O--]", "[OH-]", "[Mn++]", "[La+3]", "[Lu+3]", "[SH-]", "[Pt+4]",
               "[Fe++]", "[W]", "[Cu+2]", "[Cr+3]", "[Tc+7]", "[Xe]", "[Tl+]", "[Zn+2]", "[F-]", "[C]", "[He]", "N#N",
               "O=O", "Cl[Ra]Cl", "[Mn+2]", "N#[N+][O-]", "II", "[Ga+3]", "[Mo+10]", "[Zn]", "[Fe]", "[Si+4]", "[Al]",
               "[B+3]", "[Co]", "[Ni+2]", "[223Ra]", "[Sm]", "[Ac]", "[Mg++]", "[Fe]", "[Zn]", "[Cu]", "[Mo+5]", "[Nd]",
               "[Ta]","[Co+3]", "O", "[Cu+4]", "[Zn+4]", "[Mo]", "[Fe+4]", "[177Lu]", "[223Ra+2]", "[Ta+]", "[Cd]",
               "[153Sm]", "[Cs]", "[201Tl+]", "[V+5]", "[Hg+2]", "[Se]", "[Ba]", "[133Xe]", "[60Co]", "[98Tc-]", "[Sn+2]",
               '[51Cr+6]', "[103Pd]", "[Sb]", "[Bk]", "[Ti]", "[3He]", "[S-2]", "[90Y]", "[Co+2]", "[32PH3]", "[111In]",
               "[In]", "[Mn+7]", "[Pt]", "[201Tl]"]

uncharger = rdMolStandardize.Uncharger()



def prepInput(input):

    if input == "--" or input == "-" or input == None:
        return "Error: no valid SMILES"

    if not search(r"([a-z]|[A-Z])", input):
        SMIin = CASDTXIDtoSMILES(input)
        if SMIin == "Error":
            return "Error: CASRN not corresponding to a structure"
        else:
            return SMIin

    # Check if DSSTOX
    elif search(r"DTXSID", input.upper()):
        SMIin = CASDTXIDtoSMILES(input)
        if SMIin == "Error":
            return "Error: DSSTOX not corresponding to a structure"
        else:
            return SMIin

    return input



def prepSMI(SMIin, defnFilename="", removeMetal = 1):

    
    #step:
    try: # to hundle [H][N]([H])([H])[Pt](Cl)(Cl)[N]([H])([H])[H]
        # 1. smiles to mol
        mol = Chem.MolFromSmiles(SMIin)

        #2. normalize + clean structure
        mol_normalize = rdMolStandardize.Normalize(mol)
        mol_normalize = rdMolStandardize.Cleanup(mol_normalize)

        #3. remove metal
        if defnFilename != "":
            remover = SaltRemover(defnFilename=defnFilename)
        else:
            remover = SaltRemover()

        if removeMetal == 1:
            mol_normalize = remover.StripMol(mol_normalize)


        #4. uncharge
        mol_neutral = uncharger.uncharge(mol_normalize)

        #5. remove H
        mol_neutral_withoutH =  Chem.RemoveHs(mol_neutral)

        #6. canonical SMILES
        smilesclean = Chem.MolToSmiles(mol_neutral_withoutH,isomericSmiles=False)
    except:
        smilesclean = ""
    
    #7. mixture and other ion - inorganic elements to remove
    # 2. SMILES remove other manual salts + fragments -> for fragment take one if exactly same compound
    lelem = smilesclean.split(".")
    # reduce double, case of several salts are included - 255
    lelem = list(set(lelem))
    try:lelem.remove("")
    except:pass

    # remove metal
    if removeMetal == 1:
        lnometal = []
        for elem in lelem:
            if is_metalorion(elem) == 0:
                lnometal.append(elem)
        lelem = lnometal


    if len(lelem) == 1:
        smilesclean = str(lelem[0])
        return smilesclean
    elif len(lelem) > 1:
        return "Mixture or fragment to check: " + smilesclean
    elif smilesclean == "":
        return "Error: SMILES empty after preparation"
    else:
        return "Error: No identified"
    

# remove metal or ion
def is_metalorion(smilesin):
    # case only atom in
    if len(smilesin) <= 3:
        return 1
    if smilesin[0] == "[" and smilesin[-1] == "]":
        if len(smilesin) < 10:
            return 1
    return 0



def CASDTXIDtoSMILES(smiIn):

    try:handle = urllib.request.urlopen('https://comptox.epa.gov/dashboard/dsstoxdb/results?search=' + str(smiIn))
    except: return "Error"
    resp = handle.read()
    resp = str(resp.decode("utf8"))
    handle.close()
    try:
        SMILES = resp.split("SMILES: ")[1][0:1000]
        SMILES = SMILES.split("\\nInChI")[0]
        if SMILES == None:
            return "Error"
        else:
            return str(SMILES)
    except:
        return "Error"
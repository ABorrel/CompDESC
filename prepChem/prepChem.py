from re import search
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from molvs import Standardizer

# frag not considered
LSMILESREMOVE=["[C-]#N", "[Al+3]", "[Gd+3]", "[Pt+2]", "[Au+3]", "[Bi+3]", "[Al]", "[Si+4]", "[Fe]", "[Zn]", "[Fe+2]",
               "[Ru+8]", "[Fe+]", "[Sr++]", "[Fe+3]", "[O--]", "[OH-]", "[Mn++]", "[La+3]", "[Lu+3]", "[SH-]", "[Pt+4]",
               "[Fe++]", "[W]", "[Cu+2]", "[Cr+3]", "[Tc+7]", "[Xe]", "[Tl+]", "[Zn+2]", "[F-]", "[C]", "[He]", "N#N",
               "O=O", "Cl[Ra]Cl", "[Mn+2]", "N#[N+][O-]", "II", "[Ga+3]", "[Mo+10]", "[Zn]", "[Fe]", "[Si+4]", "[Al]",
               "[B+3]", "[Co]"]

def prepInput(input):

    # CAS-ID
    if not search(r"([a-z]|[A-Z])", input):
        SMIin = CASDTXIDtoSMILES()
        if SMIin == "ERROR":
            return "ERROR: CASRN not corresponding to a structure"

    # Check if DSSTOX
    elif search(r"DTXSID", input.upper()):
        SMIin = CASDTXIDtoSMILES()
        if SMIin == "ERROR":
            return "ERROR: DSSTOX not corresponding to a structure"

    else:
        return input



def prepSMI(SMIin):

    mol = Chem.MolFromSmiles(SMIin)
    s = Standardizer()

    try:
        molstandardized = s.standardize(mol)
        smilestandadized = Chem.MolToSmiles(molstandardized)
    except:
        return "ERROR: Standardization Fail"


    # remove salt
    # 1.default
    remover = SaltRemover()
    molclean = remover(molstandardized)
    smilesclean = Chem.MolToSmiles(molclean)

    # 2. SMILES remove other manual salts + fragments -> for fragment take one if exactly same compound
    lelem = smilesclean.split(".")
    if len(lelem) > 1:
        # reduce double, case of several salts are included - 255
        lelem = list(set(lelem))
        for smilesdel in LSMILESREMOVE:
            if smilesdel in lelem:
                lelem.remove(smilesdel)
        try:
            lelem.remove("")  # case of bad smile
        except:
            pass
        if len(lelem) == 1:
            smilesclean = str(lelem[0])
        else:

            # 3. Fragments - mixture
            # Case of fragment -> stock in log file, check after to control
            return "ERROR: Mixture or fragment ot check: " + smilesclean

        if smilesclean == "":
            return "ERROR: SMILES empty after preparation"

    return smilesclean



 def CASDTXIDtoSMILES(smiIn):

    with requests.Session() as s:
        resp = s.get('https://comptox.epa.gov/dashboard/dsstoxdb/results?search=' + str(smiIn))
        #result = str(resp)
        resp = str(resp.text)
        try:
            SMILES = resp.split("SMILES: ")[1][0:1000]
            SMILES = SMILES.split("\\nInChI")[0]
            return SMILES
        except:
            return "ERROR"
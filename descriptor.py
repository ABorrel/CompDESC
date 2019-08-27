#2D descriptors
from Desc1D2D import constitution
from Desc1D2D import molproperty
from Desc1D2D import topology
from Desc1D2D import connectivity
from Desc1D2D import kappa

#from .constitution import GetConstitutional
#from .molproperty import GetMolecularProperty
#from .topology import GetTopology
#from .connectivity import GetConnectivity
#from .kappa import GetKappa
#from .bcut import GetBurden
#from .basak import Getbasak
#from .estate import GetEstate
#from .moreaubroto import GetMoreauBrotoAuto
#from .moran import GetMoranAuto
#from .geary import GetGearyAuto
#from .charge import GetCharge
#from .moe import GetMOE

#3D descriptors
#from .geo3D import GetGeo3D
#from .cpsa3D import GetCPSA3D
#from .rdf3D import GetRdf3D
#from .morse3D import GetMorse3D
#from .whim3D import GetWhim3D


#from rdkit.Chem import Descriptors
from copy import deepcopy
from rdkit import Chem


import subprocess
from os import path, remove

# list of descriptor 2D
lCONSTITUTION = ['Weight','AWeight','nhyd','nhal','nhet','nhev','ncof','ncocl','ncobr','ncoi','ncarb','nphos','nsulph',
                 'noxy','nnitro','nring','nrot','ndonr','naccr','nsb','ndb','naro','ntb','nta','PC1','PC2','PC3','PC4',
                 'PC5','PC6']
LCOMPO = ["nheavy"]
LMOLPROP = ['LogP', 'LogP2', 'MR', 'TPSA', 'Hy', 'UI']
LTOPOLOGY = ['W','AW','J','Tigdi','Xu','GMTI','Pol','DZ','Ipc','BertzCT','Thara','Tsch','ZM1','ZM2','MZM1','MZM2',
             'Qindex','Platt','diametert','radiust','petitjeant','Sito','Hato','Geto','Arto']
LCONNECTIVITY = ['Chi0','Chi1','mChi1','Chi2','Chi3','Chi4','Chi5','Chi6','Chi7','Chi8','Chi9','Chi10','Chi3c','Chi4c',
                 'Chi4pc','Chi3ch','Chi4ch','Chi5ch','Chi6ch','knotp','Chiv0','Chiv1','Chiv2','Chiv3','Chiv4','Chiv5',
                 'Chiv6','Chiv7','Chiv8','Chiv9','Chiv10','dchi0','dchi1','dchi2','dchi3','dchi4','Chiv3c','Chiv4c',
                 'Chiv4pc','Chiv3ch','Chiv4ch','Chiv5ch','Chiv6ch','knotpv']
LKAPA = ['kappa1', 'kappa2', 'kappa3', 'kappam1', 'kappam2', 'kappam3', 'phi']
LBUCUT =["bcutp16","bcutp15","bcutp14","bcutp13","bcutp12","bcutp11","bcutp10",
        "bcutp9","bcutp8","bcutp7","bcutp6","bcutp5","bcutp4","bcutp3",
        "bcutp2","bcutp1"]
LBASAK = ['CIC0','CIC1','CIC2','CIC3','CIC4','CIC5','CIC6','SIC0','SIC1','SIC2','SIC3','SIC4','SIC5','SIC6','IC0','IC1',
          'IC2','IC3','IC4','IC5','IC6']
LESTATE = ['Smax38', 'Smax39', 'Smax34', 'Smax35', 'Smax36', 'S43', 'Smax30', 'Smax31', 'Smax32', 'Smax33', 'S57',
           'S56', 'S55', 'S54', 'S53', 'S52', 'S51', 'S50', 'Smin49', 'S59', 'S58', 'Smin69', 'Smin68', 'Smin27',
           'Sfinger30', 'Sfinger31', 'Sfinger32', 'Sfinger33', 'Sfinger34', 'Sfinger35', 'Sfinger36', 'Sfinger37',
           'Sfinger38', 'Sfinger39', 'Smax2', 'Smax3', 'Smax4', 'Smax5', 'Smax6', 'Smax7', 'Smin77', 'Smax29', 'Smax37',
           'Smax23', 'Smax22', 'Smax21', 'Smax20', 'Smax27', 'Smax26', 'Smax25', 'Smax24', 'S44', 'S45', 'S46', 'S47',
           'S40', 'S41', 'S42', 'S17', 'Smin44', 'S48', 'S49', 'Smin8', 'Smin29', 'Smin28', 'Sfinger45', 'Sfinger44',
           'Sfinger47', 'Sfinger46', 'Sfinger41', 'Sfinger40', 'Sfinger43', 'Sfinger42', 'Smax47', 'Smin73', 'Smin70',
           'Smin71', 'Sfinger49', 'Sfinger48', 'Smin74', 'Smin75', 'Smin67', 'Smin6', 'Smin9', 'Smin7', 'Smin47',
           'Smax41', 'S79', 'S78', 'Smin19', 'Smax58', 'Smax59', 'S71', 'S70', 'S73', 'S72', 'S75', 'S74', 'S77',
           'S76', 'Smax73', 'Smin78', 'Sfinger56', 'Sfinger57', 'Sfinger54', 'Sfinger55', 'Sfinger52', 'Sfinger53',
           'Sfinger50', 'Sfinger51', 'Smin61', 'Smin60', 'Smin63', 'Smin62', 'Smin65', 'Smin64', 'Sfinger58',
           'Sfinger59', 'Smin48', 'Smin42', 'Smin76', 'Smin41', 'Smin72', 'Smax40', 'Smin40', 'Smax49', 'Smax48',
           'S68', 'S69', 'S66', 'S67', 'S64', 'S65', 'S62', 'S63', 'S60', 'S61', 'Smin54', 'Smax52', 'Sfinger69',
           'Sfinger68', 'Smin50', 'Smin51', 'Smin52', 'Smin53', 'Sfinger63', 'Sfinger62', 'Sfinger61', 'Sfinger60',
           'Sfinger67', 'S10', 'Sfinger65', 'Sfinger64', 'S13', 'S12', 'Sfinger76', 'Smin56', 'S9', 'S8', 'S3', 'S2',
           'S1', 'Smin55', 'S7', 'S6', 'S5', 'S4', 'Smax78', 'Smax45', 'Smax11', 'Sfinger72', 'Smin66', 'Smax44',
           'Smax70', 'Smax71', 'Smax72', 'S14', 'Smax74', 'Smax75', 'Smax76', 'Smax77', 'Smin43', 'Smax8', 'S19',
           'S18', 'Sfinger78', 'Sfinger79', 'Smin45', 'Smax9', 'Sfinger74', 'Sfinger75', 'S11', 'Sfinger77',
           'Sfinger70', 'Sfinger71', 'S15', 'Sfinger73', 'Smax43', 'Smin16', 'Smax42', 'Smax53', 'Smax66', 'Smax65',
           'Smax64', 'Smax63', 'Smax62', 'Smax61', 'Smax60', 'Smin26', 'Smax69', 'Smax68', 'Smax0', 'Smin57', 'Smax1',
           'Smin17', 'Smin36', 'Smin37', 'Smin34', 'Smin35', 'Smin32', 'Smin33', 'Smin30', 'Smin31', 'Smax67', 'Smin46',
           'Smax51', 'Smin38', 'Smin39', 'Smax12', 'Smax13', 'Smax10', 'S16', 'Smax16', 'Smax17', 'Smax14', 'Smax15',
           'Smin20', 'Smax18', 'Smax19', 'Sfinger66', 'Smax56', 'Smax28', 'Smax57', 'Smax54', 'Smin58', 'Smax55', 'S39',
           'S38', 'Smax46', 'S35', 'S34', 'S37', 'S36', 'S31', 'S30', 'S33', 'S32', 'Smin25', 'Smin24', 'Sfinger18',
           'Sfinger19', 'Smin21', 'Smax50', 'Smin23', 'Smin22', 'Sfinger12', 'Sfinger13', 'Sfinger10', 'Sfinger11',
           'Sfinger16', 'Sfinger17', 'Sfinger14', 'Sfinger15', 'Sfinger8', 'Sfinger9', 'Smin4', 'Smin5', 'Smin2',
           'Smin3', 'Smin0', 'Smin1', 'Sfinger1', 'Sfinger2', 'Sfinger3', 'Sfinger4', 'Sfinger5', 'Sfinger6',
           'Sfinger7', 'S22', 'S23', 'S20', 'S21', 'S26', 'S27', 'S24', 'S25', 'Smin59', 'S28', 'S29', 'Smin18',
           'Smin10', 'Smin11', 'Smin12', 'Smin13', 'Smin14', 'Smin15', 'Sfinger29', 'Sfinger28', 'Sfinger27',
           'Sfinger26', 'Sfinger25', 'Sfinger24', 'Sfinger23', 'Sfinger22', 'Sfinger21', 'Sfinger20']
LMOREAUBROTO = ['ATSe1', 'ATSe2', 'ATSe3', 'ATSe4', 'ATSe5', 'ATSe6', 'ATSe7', 'ATSe8', 'ATSp8', 'ATSp3', 'ATSv8',
                'ATSp1', 'ATSp7', 'ATSp6', 'ATSp5', 'ATSp4', 'ATSv1', 'ATSp2', 'ATSv3', 'ATSv2', 'ATSv5', 'ATSv4',
                'ATSv7', 'ATSv6', 'ATSm8', 'ATSm1', 'ATSm2', 'ATSm3', 'ATSm4', 'ATSm5', 'ATSm6', 'ATSm7']
LMORAN = ['MATSv8', 'MATSp4', 'MATSp8', 'MATSv1', 'MATSp6', 'MATSv3', 'MATSv2', 'MATSv5', 'MATSv4', 'MATSv7', 'MATSv6',
          'MATSm8', 'MATSp1', 'MATSm4', 'MATSm5', 'MATSm6', 'MATSm7', 'MATSm1', 'MATSm2', 'MATSm3', 'MATSe4', 'MATSe5',
          'MATSe6', 'MATSe7', 'MATSe1', 'MATSe2', 'MATSe3', 'MATSe8', 'MATSp3', 'MATSp7', 'MATSp5', 'MATSp2']
LGEARY = ['GATSp8', 'GATSv3', 'GATSv2', 'GATSv1', 'GATSp6', 'GATSv7', 'GATSv6', 'GATSv5', 'GATSv4', 'GATSe2', 'GATSe3',
          'GATSv8', 'GATSe6', 'GATSe7', 'GATSe4', 'GATSe5', 'GATSp5', 'GATSp4', 'GATSp7', 'GATSe1', 'GATSp1', 'GATSp3',
          'GATSp2', 'GATSe8', 'GATSm2', 'GATSm3', 'GATSm1', 'GATSm6', 'GATSm7', 'GATSm4', 'GATSm5', 'GATSm8']
LCHARGE = ['SPP','LDI','Rnc','Rpc','Mac','Tac','Mnc','Tnc','Mpc','Tpc','Qass','QOss','QNss','QCss','QHss','Qmin','Qmax',
           'QOmin','QNmin','QCmin','QHmin','QOmax','QNmax','QCmax','QHmax']
LMOE = ['EstateVSA8', 'EstateVSA9', 'EstateVSA4', 'EstateVSA5', 'EstateVSA6', 'EstateVSA7', 'EstateVSA0', 'EstateVSA1',
        'EstateVSA2', 'EstateVSA3', 'PEOEVSA13', 'PEOEVSA12', 'PEOEVSA11', 'PEOEVSA10', 'MTPSA', 'VSAEstate0',
        'VSAEstate1', 'VSAEstate2', 'VSAEstate3', 'VSAEstate4', 'VSAEstate5', 'VSAEstate6', 'VSAEstate7', 'VSAEstate8',
        'LabuteASA', 'PEOEVSA3', 'PEOEVSA2', 'PEOEVSA1', 'PEOEVSA0', 'PEOEVSA7', 'PEOEVSA6', 'PEOEVSA5', 'PEOEVSA4',
        'MRVSA5', 'MRVSA4', 'PEOEVSA9', 'PEOEVSA8', 'MRVSA1', 'MRVSA0', 'MRVSA3', 'MRVSA2', 'MRVSA9', 'slogPVSA10',
        'slogPVSA11', 'MRVSA8', 'MRVSA7', 'MRVSA6', 'EstateVSA10', 'slogPVSA2', 'slogPVSA3', 'slogPVSA0', 'slogPVSA1',
        'slogPVSA6', 'slogPVSA7', 'slogPVSA4', 'slogPVSA5', 'slogPVSA8', 'slogPVSA9', 'VSAEstate9']#, 'VSAEstate10']

L3D = ['RDFC6', 'MoRSEN11', 'RDFU8', 'RDFU9', 'RDFU2', 'RDFU3', 'MoRSEN5', 'RDFU1', 'RDFU6', 'RDFU7', 'RDFU4', 'RDFU5',
       'Harary3D', 'P2u', 'MoRSEM6', 'MoRSEM7', 'MoRSEM4', 'MoRSEM5', 'MoRSEM2', 'MoRSEM3', 'MoRSEE30', 'MoRSEM1',
       'MoRSEN4', 'MoRSEM8', 'MoRSEM9', 'MoRSEU10', 'MoRSEU11', 'MoRSEU12', 'MoRSEU13', 'MoRSEU14', 'MoRSEU15',
       'MoRSEU16', 'MoRSEU17', 'MoRSEU18', 'MoRSEU19', 'FPSA3', 'FPSA2', 'FPSA1', 'GeDi', 'MoRSEV19', 'MoRSEN10',
       'MoRSEV13', 'SPAN', 'MoRSEV11', 'MoRSEV10', 'MoRSEV17', 'MoRSEV16', 'MoRSEV15', 'MoRSEN16', 'RDFM14', 'RDFM15',
       'RDFM16', 'RDFE9', 'RDFM10', 'RDFM11', 'RDFM12', 'RASA', 'RDFE2', 'RDFE3', 'MoRSEC25', 'MoRSEC24', 'RDFM18',
       'RDFM19', 'grav', 'RDFE5', 'WNSA1', 'WNSA2', 'WNSA3', 'L2p', 'RDFP15', 'RDFP14', 'RDFP17', 'RDFP16', 'RDFP11',
       'RDFP10', 'RDFP13', 'RDFP12', 'MoRSEP30', 'RDFP19', 'RDFP18', 'E2p', 'Dm', 'P3e', 'MoRSEM18', 'MoRSEM19',
       'Petitj3D', 'MoRSEM10', 'MoRSEM11', 'MoRSEM12', 'MoRSEM13', 'MoRSEM14', 'MoRSEM15', 'MoRSEM16', 'MoRSEM17',
       'RDFC27', 'RDFC26', 'RDFC25', 'RDFC24', 'RDFC23', 'RDFC22', 'RDFC21', 'RDFC20', 'MoRSEC28', 'RDFU30', 'RDFC29',
       'RDFC28', 'MoRSEU30', 'L1u', 'L1v', 'L2v', 'L1p', 'RDFP5', 'RDFP4', 'RDFP7', 'RDFP6', 'RDFP1', 'RDFP3', 'RDFP2',
       'L1e', 'RDFP9', 'RDFP8', 'MoRSEP5', 'P1e', 'MoRSEP4', 'PSA', 'MoRSEP7', 'P1p', 'MoRSEP6', 'RDFE18', 'RDFE19',
       'P1v', 'RDFE14', 'RDFE15', 'RDFE16', 'RDFE17', 'PPSA1', 'RDFE11', 'RDFE12', 'PPSA2', 'MoRSEP11', 'MoRSEP10',
       'MoRSEP13', 'MoRSEP12', 'RPCS', 'MoRSEP14', 'MoRSEN9', 'MoRSEN8', 'DPSA1', 'MoRSEC30', 'DPSA3', 'DPSA2',
       'MoRSEN3', 'MoRSEN2', 'MoRSEN1', 'RDFP30', 'E2e', 'MoRSEN17', 'L3e', 'TASA', 'RDFC19', 'MoRSEV14', 'MoRSEM30',
       'MoRSEP8', 'L3v', 'RDFC16', 'L3u', 'RDFV30', 'L3p', 'RDFC14', 'W3DH', 'RDFC15', 'MoRSEC23', 'MoRSEN15',
       'MoRSEP16', 'RPSA', 'P3m', 'MEcc', 'MoRSEC22', 'MoRSEN14', 'MoRSEP1', 'MoRSEN23', 'P3p', 'P3v', 'MoRSEP19',
       'P3u', 'RDFV7', 'RDFC18', 'RDFV6', 'FNSA1', 'RDFC17', 'FNSA3', 'FNSA2', 'RDFC12', 'RDFC13', 'RDFC10', 'RDFC11',
       'P2p', 'RDFV4', 'MoRSEP22', 'RDFV3', 'MoRSEP18', 'RDFV2', 'RDFU21', 'RDFU20', 'RDFU23', 'RDFU22', 'RDFU25',
       'RDFM4', 'RDFU27', 'RDFU26', 'RDFU29', 'RDFU28', 'MoRSEN28', 'MoRSEN29', 'RDFV19', 'RDFV18', 'WPSA2', 'RDFV16',
       'RDFV15', 'RDFV14', 'RDFV13', 'RDFV12', 'RDFV11', 'RDFV10', 'MoRSEN26', 'MoRSEP23', 'MoRSEN27', 'MoRSEN24',
       'MoRSEP25', 'MoRSEN25', 'MoRSEP26', 'MoRSEE8', 'MoRSEE9', 'MoRSEE6', 'MoRSEN22', 'MoRSEE4', 'MoRSEP27',
       'MoRSEE2', 'MoRSEE3', 'RDFU24', 'MoRSEN20', 'MoRSEC9', 'ASPAN', 'RDFE10', 'MoRSEN21', 'Te', 'Vm', 'Vp',
       'MoRSEV18', 'PPSA3', 'Vv', 'RDFE13', 'E2u', 'RDFC30', 'E2v', 'P1m', 'MoRSEV12', 'MoRSEP15', 'MoRSEP17',
       'MoRSEU8', 'MoRSEU9', 'MoRSEU6', 'MoRSEU7', 'MoRSEU4', 'MoRSEU5', 'MoRSEU2', 'MoRSEU3', 'MoRSEU1', 'RDFM29',
       'RDFM28', 'MoRSEN6', 'RDFM21', 'RDFM20', 'RDFM23', 'RDFM22', 'RDFM25', 'RDFM24', 'RDFM27', 'RDFM26', 'RDFM2',
       'RDFM3', 'RDFV5', 'RDFM1', 'RDFM6', 'RDFM7', 'RDFV1', 'RDFM5', 'MoRSEP20', 'MoRSEP21', 'RDFM8', 'RDFM9',
       'MoRSEP24', 'RDFE8', 'RDFV9', 'RDFV8', 'MoRSEC8', 'RDFV17', 'RDFM17', 'WPSA3', 'AGDD', 'MoRSEC1', 'MoRSEC2',
       'MoRSEC3', 'MoRSEC4', 'MoRSEC5', 'MoRSEC6', 'MoRSEC7', 'MoRSEE29', 'MoRSEE28', 'WPSA1', 'MoRSEC29', 'MoRSEE21',
       'MoRSEE20', 'MoRSEE23', 'RDFM13', 'MoRSEE25', 'MoRSEE24', 'MoRSEE27', 'MoRSEE26', 'MoRSEC27', 'MoRSEC26', 'Ae',
       'RDFE1', 'RDFE6', 'RDFE7', 'RDFE4', 'MoRSEV28', 'MoRSEV29', 'MoRSEV26', 'MoRSEV27', 'MoRSEV24', 'MoRSEC20',
       'MoRSEV22', 'MoRSEV23', 'MoRSEV20', 'MoRSEV21', 'MoRSEC12', 'MoRSEC13', 'MoRSEC10', 'MoRSEC11', 'MoRSEC16',
       'MoRSEC17', 'MoRSEC14', 'MoRSEC15', 'P2m', 'MoRSEC18', 'MoRSEC19', 'RDFE30', 'RDFE21', 'RDFE20', 'RDFE23',
       'RDFE22', 'RDFE25', 'RDFE24', 'RDFE27', 'RDFE26', 'RDFE29', 'RDFE28', 'RDFP20', 'RDFP21', 'RDFP22', 'RDFP23',
       'RDFP24', 'RDFP25', 'RDFP26', 'RDFP27', 'RDFP28', 'RDFP29', 'MoRSEV6', 'MoRSEN13', 'MoRSEV30', 'Dv', 'RDFV26',
       'RDFV27', 'RDFV24', 'RDFV25', 'RDFV22', 'RDFV23', 'RDFV20', 'RDFV21', 'L2e', 'MoRSEV5', 'RDFV28', 'RDFV29',
       'MoRSEP3', 'P1u', 'rygr', 'Ve', 'MoRSEE7', 'MoRSEV4', 'MoRSEP2', 'FrTATP', 'MoRSEE5', 'P2v', 'ASA', 'MoRSEC21',
       'MoRSEV3', 'MoRSEE1', 'E3v', 'E3u', 'Ke', 'E3p', 'Km', 'MoRSEV7', 'E3e', 'Kp', 'Kv', 'Ku', 'MoRSEV2', 'RDFC8',
       'E3m', 'RDFC9', 'MoRSEM21', 'MoRSEM20', 'MoRSEM23', 'MoRSEM22', 'MoRSEM25', 'MoRSEM24', 'MoRSEM27', 'MoRSEM26',
       'MoRSEM29', 'MoRSEM28', 'MoRSEN19', 'MoRSEN18', 'L2m', 'MoRSEV9', 'MoRSEV8', 'MoRSEU29', 'MoRSEU28', 'L2u',
       'MoRSEV1', 'MoRSEN12', 'RDFC5', 'MoRSEU21', 'MoRSEU20', 'MoRSEU23', 'MoRSEU22', 'MoRSEU25', 'MoRSEU24',
       'MoRSEU27', 'MoRSEU26', 'RDFC7', 'MoRSEV25', 'Tv', 'Am', 'SEig', 'Tu', 'Tp', 'RDFC1', 'Tm', 'RDFC2', 'PNSA3',
       'PNSA2', 'PNSA1', 'RDFC3', 'MSA', 'MoRSEE22', 'L1m', 'De', 'MoRSEE10', 'MoRSEE11', 'MoRSEE12', 'MoRSEE13',
       'MoRSEP9', 'MoRSEE15', 'MoRSEE16', 'MoRSEE17', 'MoRSEE18', 'MoRSEE19', 'Du', 'Dp', 'Vu', 'P2e', 'E1p', 'E1u',
       'E1v', 'E1m', 'RNCS', 'MoRSEP28', 'MoRSEN30', 'E1e', 'MoRSEN7', 'W3D', 'RDFU18', 'RDFU19', 'RDFU14', 'RDFU15',
       'RDFU16', 'RDFU17', 'RDFU10', 'RDFU11', 'RDFU12', 'RDFU13', 'Ap', 'Au', 'RDFC4', 'MoRSEP29', 'Av', 'L3m',
       'RDFM30', 'MoRSEE14', 'E2m']



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


def getLdesc (typeDesc):

    lout = []
    if typeDesc == "1D2D":
        # listdesc
        lout = LCOMPO + LMOLPROP + LTOPOLOGY + LCONNECTIVITY + LKAPA + LBUCUT + LBASAK + LESTATE + LMOREAUBROTO + LMORAN + LGEARY + LCHARGE + LMOE

    elif typeDesc == "3D":
        lout = L3D

    return lout




def computePNG(SMILES, inchikey, prin, prout):

    pSMILES = prin + inchikey + ".smi"
    pPNG = prout + inchikey + ".png"
    if path.exists(pPNG):
        return
    else:
        fSMI = open(pSMILES, "w")
        fSMI.write(str(SMILES))
        fSMI.close()
        cmd = "molconvert \"png:w500,Q100,#00000000\" " + pSMILES + " -o " + pPNG
        subprocess.Popen(cmd, shell=True)


class Descriptor:

    # mol have to be clean before
    def __init__(self, SMICLEAN, pdesc):
        self.smi = SMICLEAN
        self.mol = Chem.MolFromSmiles(SMICLEAN)
        self.err = 0
        self.pdesc = pdesc

    def computeAll2D(self):

        if path.exists(self.pdesc + "_2D.txt"):
            if path.getsize(self.pdesc + "_2D.txt") > 100:
                ddesc = loadMatrixToDict(self.pdesc + "_2D.txt")
                self.all2D = ddesc
                return
            else:
                remove(self.pdesc + "_2D.txt")

        #self.consti = constitution.GetConstitutional(self.mol)
        #self.molprop = molproperty.GetMolecularProperty(self.mol)
        #self.topology = topology.GetTopology(self.mol)
        #self.connectivity = connectivity.GetConnectivity(self.mol)
        self.kappa = kappa.GetKappa(self.mol)
        return



        self.burden = GetBurden(self.mol)
        self.basakD = Getbasak(self.mol)
        self.estate = GetEstate(self.mol)
        self.moreauBurto = GetMoreauBrotoAuto(self.mol)
        self.autcormoran = GetMoranAuto(self.mol)
        self.gearycor = GetGearyAuto(self.mol)
        self.charges = GetCharge(self.mol)
        self.MOE = GetMOE(self.mol)

        # combine 2D
        self.all2D = {}
        self.all2D.update(deepcopy(self.consti))
        self.all2D.update(deepcopy(self.compo))
        self.all2D.update(deepcopy(self.molprop))
        self.all2D.update(deepcopy(self.topology))
        self.all2D.update(deepcopy(self.connectivity))
        self.all2D.update(deepcopy(self.kappa))
        self.all2D.update(deepcopy(self.burden))
        self.all2D.update(deepcopy(self.basakD))
        self.all2D.update(deepcopy(self.estate))
        self.all2D.update(deepcopy(self.moreauBurto))
        self.all2D.update(deepcopy(self.autcormoran))
        self.all2D.update(deepcopy(self.gearycor))
        self.all2D.update(deepcopy(self.charges))
        self.all2D.update(deepcopy(self.MOE))




    def computeAll3D(self,lcoordinates):

        if path.exists(self.pdesc + "_3D.txt"):
            if path.getsize(self.pdesc + "_3D.txt") > 100:
                self.all3D = loadMatrixToDict(self.pdesc + "_3D.txt")
                return
            else:
                remove(self.pdesc + "_3D.txt")


        self.geo3D = GetGeo3D(lcoordinates)
        self.CPSA3D = GetCPSA3D(lcoordinates)
        self.rdf3D = GetRdf3D(lcoordinates)
        self.morse3D = GetMorse3D(lcoordinates)
        self.whim3D = GetWhim3D(lcoordinates)

        # combine 3D
        self.all3D = {}
        self.all3D.update(deepcopy(self.geo3D))
        self.all3D.update(deepcopy(self.CPSA3D))
        self.all3D.update(deepcopy(self.rdf3D))
        self.all3D.update(deepcopy(self.morse3D))
        self.all3D.update(deepcopy(self.whim3D))


    def writeMatrix(self, typedesc):
        if typedesc == "2D":
            if "all2D" in self.__dict__:
                filin = open(self.pdesc + "_2D.txt", "w")
                filin.write("%s\n"%("\t".join(self.all2D.keys())))
                filin.write("%s\n"%("\t".join([str(self.all2D[k]) for k in self.all2D.keys()])))
                filin.close()

        if typedesc == "3D":
            if "all3D" in self.__dict__:
                filin = open(self.pdesc + "_3D.txt", "w")
                filin.write("%s\n"%("\t".join(self.all3D.keys())))
                filin.write("%s\n"%("\t".join([str(self.all3D[k]) for k in self.all3D.keys()])))
                filin.close()



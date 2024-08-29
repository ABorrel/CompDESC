import unittest
import CompDesc
from shutil import rmtree
from os import system


class TestChemical(unittest.TestCase):
        
    def test_PrepChem(self):
        cChem = CompDesc.CompDesc("O=S(=O)(O)c1ccc([Hg])cc1", "./tests/")
        cChem.prepChem()
        smi = cChem.smi
        self.assertEqual(smi, "O=S(=O)(O)c1ccc([Hg])cc1")

    def test_mixturePrep(self):
        cChem = CompDesc.CompDesc("ClC1=CC2=C(C=C1Cl)C1=C(O2)C(Cl)=C(Cl)C(Cl)=C1Cl.ClC1=CC2=C(C3=C(O2)C(Cl)=C(Cl)C(Cl)=C3)C(Cl)=C1Cl", "./tests/")
        cChem.prepChem()
        smi = cChem.smi
        self.assertEqual(cChem.isfrag, 1)

    def test_compute1D2Ddesc(self):
        cChem = CompDesc.CompDesc("N=C(O)[C@@H](N)CS", "./tests/")
        cChem.prepChem()
        cChem.computeAll2D()
        self.assertEqual(cChem.err, 0)

    def test_knimeConvert(self):
        cChem = CompDesc.CompDesc("N=C(O)[C@@H](N)CS", "./tests/")
        cChem.prepChem()
        cChem.computeAll2D()
        cChem.convertDesc2DtoKnimeDesc()
        err = 0
        try: test = cChem.all2D["AMW"]
        except:err = 1
        self.assertEqual(err, 0)

    def test_generate3D(self):
        cChem = CompDesc.CompDesc("N=C(O)[C@@H](N)CS", "./tests/")
        cChem.prepChem()
        cChem.set3DChemical()
        self.assertEqual(cChem.err, 0)
        rmtree("./tests/MOLCLEAN")
        rmtree("./tests/SDF3D")

    def test_compute3Ddesc(self):
        cChem = CompDesc.CompDesc("N=C(O)[C@@H](N)C", "./tests/")
        cChem.prepChem()
        cChem.set3DChemical()
        cChem.computeAll3D()
        cChem.writeMatrix("3D")
        self.assertEqual(cChem.err, 0)
        rmtree("./tests/MOLCLEAN")
        rmtree("./tests/SDF3D")
        rmtree("./tests/3D")
    
    #def test_computeOPERA(self):
    #    cChem = CompDesc.CompDesc("N=C(O)[C@@H](N)CS", "./tests/")
    #    cChem.prepChem()
    #    cChem.computeOPERAFromChem(update=1)
    #    self.assertEqual(cChem.err, 0)
    #    system("rm -rf ./tests/PADEL*")
    #    system("rm -rf ./tests/CDK*")
    #    system("rm -rf ./tests/cdk_desc*")
    #    rmtree("./tests/OPERA")

    #def test_computeOPERAServer(self):
    #    cChem = CompDesc.CompDesc("N=C(O)[C@@H](N)CS", "./tests/")
    #    cChem.prepChem()
    #    cChem.computePADEL2DFPandCDK()
    #    cChem.computeOperaDesc()
    #    self.assertEqual(cChem.err, 0)
    #    system("rm -rf ./tests/PADEL*")
    #    system("rm -rf ./tests/cdk_desc*")
    #    rmtree("./tests/OPERA")
    
    
    def test_FP(self):

        cChem = CompDesc.CompDesc("N=C(O)[C@@H](N)CS", "./tests/")
        cChem.prepChem()
        cChem.computeFP("All")
        
        # test comparison
        cCchem2 = CompDesc.CompDesc("CCCO", "./tests/")
        cCchem2.prepChem()
        cCchem2.computeFP("All")

        l_dist = ["Tanimoto", "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "Kulczynski", "McConnaughey", "Asymmetric", "BraunBlanquet", "AllBit"]
        l_fp = ['Mol', 'pairs', 'MACCS', 'Torsion', 'Morgan']
        for fp in l_fp:
            for dist in l_dist:
                print(fp, dist, cChem.computeSimilarityFP(cCchem2, fp, dist)) 
        # test all combination
        self.assertEqual(cChem.err, 0)


if __name__ == '__main__':
    unittest.main()

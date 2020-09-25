import unittest
import Chemical
from shutil import rmtree


class TestChemical(unittest.TestCase):

    def test_PrepChem(self):
        cChem = Chemical.Chemical("O=S(=O)(O)c1ccc([Hg])cc1", "./tests/")
        cChem.prepChem()
        smi = cChem.smi
        self.assertEqual(smi, "O=S(=O)(O)c1ccc([Hg])cc1")

    def test_compute1D2Ddesc(self):
        cChem = Chemical.Chemical("N=C(O)[C@@H](N)CS", "./tests/")
        cChem.prepChem()
        cChem.computeAll2D()
        self.assertEqual(cChem.err, 0)

    def test_generate3D(self):
        cChem = Chemical.Chemical("N=C(O)[C@@H](N)CS", "./tests/")
        cChem.prepChem()
        cChem.set3DChemical()
        self.assertEqual(cChem.err, 0)
        rmtree("./tests/MOLCLEAN")
        rmtree("./tests/SDF3D")

    def test_compute3Ddesc(self):
        cChem = Chemical.Chemical("CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CCCCN)C(=O)N[C@H](C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CC(C)C)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCCNC(=N)N)C(=O)O)C(C)C", "./tests/")
        cChem.prepChem()
        cChem.set3DChemical()
        cChem.computeAll3D()
        self.assertEqual(cChem.err, 0)
        rmtree("./tests/MOLCLEAN")
        rmtree("./tests/SDF3D")

    
if __name__ == '__main__':
    unittest.main()

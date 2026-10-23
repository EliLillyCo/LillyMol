import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))

import lillymol_depict as depict
import lillymol


class TestDepict(unittest.TestCase):

    def test_generate_in_place_with_keywords(self):
        mol = lillymol.MolFromSmiles("CCOc1ccccc1")
        result = depict.generate_2d_coordinates(
            mol, precision=depict.Precision.BEST, bond_length=2.0)

        self.assertIsInstance(result, depict.Coords2DResult)
        self.assertTrue(all(math.isfinite(value)
                            for value in mol.get_coordinates()))
        self.assertTrue(all(mol.z(i) == 0.0 for i in range(mol.natoms())))
        lengths = [mol.bond_length(bond.a1(), bond.a2())
                   for bond in mol.bonds()]
        self.assertAlmostEqual(sum(lengths) / len(lengths), 2.0, places=4)

    def test_options_and_copy_leave_input_unchanged(self):
        mol = lillymol.MolFromSmiles("CCCO")
        for i in range(mol.natoms()):
            mol.setxyz(i, float(i), float(i + 1), 3.0)
        original = mol.get_coordinates()

        options = depict.Coords2DOptions(centre=False)
        depicted = depict.generate_2d_coordinates_copy(mol, options)

        self.assertIsNotNone(depicted)
        copy, result = depicted
        self.assertIsInstance(result, depict.Coords2DResult)
        self.assertEqual(mol.get_coordinates(), original)
        self.assertNotEqual(copy.get_coordinates(), original)
        self.assertTrue(all(copy.z(i) == 0.0 for i in range(copy.natoms())))

    def test_empty_molecule_returns_none_without_change(self):
        mol = lillymol.Molecule()
        self.assertIsNone(depict.generate_2d_coordinates(mol))
        self.assertEqual(mol.natoms(), 0)

    def test_assign_wedge_bonds(self):
        mol = lillymol.MolFromSmiles("C[C@H](N)C(=O)O")
        self.assertIsNotNone(depict.generate_2d_coordinates(mol))

        result = depict.assign_wedge_bonds(mol)

        self.assertEqual(result.chiral_centres, 1)
        self.assertEqual(result.wedged, 1)
        self.assertEqual(result.unresolved, 0)


if __name__ == "__main__":
    unittest.main()

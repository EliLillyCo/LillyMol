import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from lillymol import Molecule
from lillymol_tools import MedchemWizard


def _reactions_file():
    candidates = []

    lillymol_home = os.environ.get("LILLYMOL_HOME")
    if lillymol_home:
        candidates.append(os.path.join(lillymol_home, "data", "MedchemWizard", "REACTIONS"))

    candidates.extend([
        os.path.join(os.getcwd(), "data", "MedchemWizard", "REACTIONS"),
        os.path.join(os.path.dirname(__file__), "..", "..", "data", "MedchemWizard", "REACTIONS"),
    ])

    test_srcdir = os.environ.get("TEST_SRCDIR")
    if test_srcdir:
        candidates.extend([
            os.path.join(test_srcdir, "_main", "data", "MedchemWizard", "REACTIONS"),
            os.path.join(test_srcdir, "data", "MedchemWizard", "REACTIONS"),
        ])

    runfiles_dir = os.environ.get("RUNFILES_DIR")
    if runfiles_dir:
        candidates.extend([
            os.path.join(runfiles_dir, "_main", "data", "MedchemWizard", "REACTIONS"),
            os.path.join(runfiles_dir, "data", "MedchemWizard", "REACTIONS"),
        ])

    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate

    return None


def _methane():
    mol = Molecule()
    assert mol.build_from_smiles("C methane")
    return mol


class TestMedchemWizard(unittest.TestCase):

    def _wizard(self):
        reactions = _reactions_file()
        if reactions is None:
            self.skipTest("Cannot locate MedchemWizard REACTIONS data")

        wizard = MedchemWizard()
        wizard.read_reactions(reactions)
        wizard.set_max_atoms(3)
        wizard.set_append_names(True)
        wizard.set_name_separator(" ")
        return wizard

    def test_initialise_from_environment(self):
        if "LILLYMOL_HOME" not in os.environ:
            self.skipTest("LILLYMOL_HOME not set")

        wizard = MedchemWizard()
        wizard.initialise_from_environment()
        self.assertGreater(wizard.number_reactions(), 0)

    def test_process_returns_products_and_does_not_change_input(self):
        wizard = self._wizard()
        mol = _methane()
        initial_smiles = mol.smiles()

        products = wizard.process(mol)

        self.assertGreater(len(products), 0)
        self.assertEqual(mol.smiles(), initial_smiles)
        self.assertTrue(all(product.natoms() <= 3 for product in products))

        stats = wizard.stats()
        self.assertEqual(stats["molecules_read"], 1)
        self.assertGreaterEqual(stats["molecules_produced"], len(products))

    def test_protected_atoms_suppress_products(self):
        wizard = self._wizard()
        wizard.add_do_not_change_smarts("[*]")

        products = wizard.process(_methane())

        self.assertEqual(len(products), 0)
        self.assertGreater(wizard.stats()["embeddings_rejected_for_changing_protected_atoms"], 0)

    def test_protected_query_no_match_can_be_ignored(self):
        wizard = self._wizard()
        wizard.add_do_not_change_smarts("[Cl]")

        with self.assertRaises(RuntimeError):
            wizard.process(_methane())

        wizard = self._wizard()
        wizard.add_do_not_change_smarts("[Cl]")
        wizard.set_ignore_do_not_change_queries_not_matching(True)
        products = wizard.process(_methane())
        self.assertGreater(len(products), 0)


if __name__ == "__main__":
    unittest.main()

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from lillymol import Molecule
from lillymol_tools import QED


def _qed_dir():
    candidates = []

    lillymol_home = os.environ.get("LILLYMOL_HOME")
    if lillymol_home:
        candidates.append(os.path.join(lillymol_home, "data", "queries", "QED"))

    candidates.extend([
        os.path.join(os.getcwd(), "data", "queries", "QED"),
        os.path.join(os.path.dirname(__file__), "..", "..", "data", "queries", "QED"),
    ])

    test_srcdir = os.environ.get("TEST_SRCDIR")
    if test_srcdir:
        candidates.extend([
            os.path.join(test_srcdir, "_main", "data", "queries", "QED"),
            os.path.join(test_srcdir, "data", "queries", "QED"),
        ])

    runfiles_dir = os.environ.get("RUNFILES_DIR")
    if runfiles_dir:
        candidates.extend([
            os.path.join(runfiles_dir, "_main", "data", "queries", "QED"),
            os.path.join(runfiles_dir, "data", "queries", "QED"),
        ])

    for candidate in candidates:
        if os.path.exists(os.path.join(candidate, "unwanted-groups.smt")):
            return candidate

    return None


class TestQED(unittest.TestCase):

    def test_qed_from_directory(self):
        mol = Molecule()
        self.assertTrue(mol.build_from_smiles("CC(=O)Oc1ccccc1C(=O)O aspirin"))
        initial_smiles = mol.smiles()

        dirname = _qed_dir()
        if dirname is None:
            self.skipTest("Cannot locate QED query data")

        qed = QED(initialise_from_environment=False)
        self.assertTrue(qed.initialise_from_directory(dirname))

        value = qed.qed(mol)
        self.assertIsNotNone(value)
        self.assertGreaterEqual(value, 0.0)
        self.assertLessEqual(value, 1.0)
        self.assertEqual(mol.smiles(), initial_smiles)


if __name__ == "__main__":
    unittest.main()

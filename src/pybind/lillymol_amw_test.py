"""Molecular weight, and what an isotopic label does to it.

An isotope in LillyMol is usually an arbitrary atom marker rather than a
statement about the nucleus. amw() therefore erases the labels before weighing -
the mass and the Hydrogen count both - so a labelled molecule weighs the same as
the unlabelled one. molecule_filter's amw feature does the same thing, and these
tests exist so that the two cannot drift apart.

amw_ignore_isotopes() deliberately differs on one point: it honours the Hydrogen
count a bracket atom declares. Both behaviours are pinned here.
"""

import os
import sys
import unittest
import warnings

sys.path.insert(0, os.path.dirname(__file__))

from lillymol import MolFromSmiles, fraction_csp3


class TestAmw(unittest.TestCase):

    # A labelled molecule must weigh the same as the unlabelled one, whatever
    # the label is. The [37C] case is the interesting one - a strict smiles
    # reading says it has no Hydrogens, which is almost never what was meant.
    SAME_WEIGHT = (
        ("COC", "[13C]OC"),
        ("COC", "[37C]OC"),
        ("COC", "[13CH3]OC"),
        ("COC", "[37CH3]OC"),
        ("CCO", "[2H]OCC"),
        ("CCO", "CC[18O]"),
        ("c1ccccc1", "[13cH]1ccccc1"),
        ("CC(=O)N", "[15NH2]C(C)=O"),
    )

    def test_label_does_not_change_amw(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for plain, labelled in self.SAME_WEIGHT:
                with self.subTest(plain=plain, labelled=labelled):
                    self.assertAlmostEqual(MolFromSmiles(plain).amw(),
                                           MolFromSmiles(labelled).amw(),
                                           places=3)

    def test_no_isotope_no_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertAlmostEqual(MolFromSmiles("COC").amw(), 46.068, places=3)

    def test_isotope_warns(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with self.assertRaises(UserWarning):
                MolFromSmiles("[37C]OC").amw()

    def test_amw_never_zero_for_isotopes(self):
        """It used to return 0.0, which reads as a plausible value."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for smiles in ("[37C]OC", "[2H]OCC", "[13cH]1ccccc1"):
                with self.subTest(smiles=smiles):
                    self.assertGreater(MolFromSmiles(smiles).amw(), 0.0)

    # amw_ignore_isotopes honours the declared Hydrogen count, so it differs
    # from amw exactly when a bracket atom declares zero Hydrogens.
    def test_amw_ignore_isotopes_honours_declared_hydrogens(self):
        for smiles, expected in (("[37C]OC", 43.045),
                                 ("[13C]OC", 43.045),
                                 ("[37CH3]OC", 46.068),
                                 ("COC", 46.068)):
            with self.subTest(smiles=smiles):
                self.assertAlmostEqual(MolFromSmiles(smiles).amw_ignore_isotopes(),
                                       expected, places=3)


class TestFractionCsp3(unittest.TestCase):

    def test_fraction_csp3(self):
        for smiles, expected in (("CCCCc1ccccc1", 0.4),   # 4 saturated of 10 carbons
                                 ("c1ccccc1", 0.0),
                                 ("CCC", 1.0),
                                 ("O", 0.0),              # no carbon at all
                                 ("CC#CC", 0.5)):
            with self.subTest(smiles=smiles):
                self.assertAlmostEqual(fraction_csp3(MolFromSmiles(smiles)),
                                       expected, places=4)


if __name__ == "__main__":
    unittest.main()

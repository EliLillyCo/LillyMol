# Tests for unique_molecules

from absl import app
from absl import logging
from absl.testing import absltest

from lillymol import *
from lillymol_tools import UniqueMolecules

class TestUniqueMolecules(absltest.TestCase):
  def test_methane(self):
    smiles = "C"
    mols = [MolFromSmiles(smi) for smi in smiles * 10]
    u = UniqueMolecules()
    self.assertTrue(u.is_unique(mols.pop()))
    for mol in mols:
      self.assertFalse(u.is_unique(mol))

  def test_chiral_different(self):
    m1 = MolFromSmiles("C(O)[C@@H](N)C")
    m2 = MolFromSmiles("C(O)[C@H](N)C")

    u = UniqueMolecules()
    self.assertTrue(u.is_unique(m1))
    self.assertTrue(u.is_unique(m2))
    
  def test_chiral_same(self):
    m1 = MolFromSmiles("C(O)[C@@H](N)C")
    m2 = MolFromSmiles("C(O)[C@H](N)C")

    u = UniqueMolecules()
    u.set_exclude_chiral_info(True)
    self.assertTrue(u.is_unique(m1))
    self.assertFalse(u.is_unique(m2))

  def test_fragments_different(self):
    m1 = MolFromSmiles("CC.C")
    m2 = MolFromSmiles("CC.O")

    u = UniqueMolecules()
    self.assertTrue(u.is_unique(m1))
    self.assertTrue(u.is_unique(m2))
    
  def test_fragments_same(self):
    m1 = MolFromSmiles("CC.C")
    m2 = MolFromSmiles("CC.O")

    u = UniqueMolecules()
    u.set_strip_to_largest_fragment(True)
    self.assertTrue(u.is_unique(m1))
    self.assertFalse(u.is_unique(m2))

  def test_etrans_not_in_effect(self):
    m1 = MolFromSmiles("Ic1cc(Cl)c(Br)c(F)c1")
    m2 = MolFromSmiles("Clc1cc(I)c(Br)c(F)c1")
    
    u = UniqueMolecules()
    self.assertTrue(u.is_unique(m1))
    self.assertTrue(u.is_unique(m2))

  def test_etrans_applied(self):
    return None
    # Temporarily disabled till we figure out the Element ptr issue.
    m1 = MolFromSmiles("Ic1cc(Cl)c(Br)c(F)c1")
    m2 = MolFromSmiles("Clc1cc(I)c(Br)c(F)c1")
    
    u = UniqueMolecules()
    u.add_element_transformation("I=Cl")
    u.add_element_transformation("Br=Cl")

    self.assertTrue(u.is_unique(m1))
    self.assertFalse(u.is_unique(m2))

  def test_isotopes_differentiated(self):
    m1 = MolFromSmiles("C")

    u = UniqueMolecules()
    self.assertTrue(u.is_unique(m1))
    self.assertFalse(u.is_unique(m1))

    m1.set_isotope(0, 1);
    self.assertTrue(u.is_unique(m1))

  def test_isotopes_same(self):
    m1 = MolFromSmiles("C")

    u = UniqueMolecules()
    u.set_ignore_isotopes(True);
    self.assertTrue(u.is_unique(m1))

    m1.set_isotope(0, 1)
    self.assertFalse(u.is_unique(m1))

if __name__ == '__main__':
  absltest.main()

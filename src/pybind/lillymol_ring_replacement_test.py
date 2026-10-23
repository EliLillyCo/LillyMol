from absl import app
from absl import logging
from absl.testing import absltest

from lillymol import *
from lillymol_tools import *

class TestRingReplacement(absltest.TestCase):
  def test_read_default(self):
    # Test that we can find replacement rings from LILLYMOL_HOME
    rr = RingReplacement()
    self.assertGreater(rr.read_replacement_rings("rings_5a5a.smi"), 0)

  def test_no_substructure_match(self):
    m = MolFromSmiles("Oc1ccc(OC)cc1")
    self.assertIsNotNone(m)

    rr = RingReplacement()
    self.assertTrue(rr.set_ring_atom_smarts("[#7]"))
    self.assertGreater(rr.read_replacement_rings("${LILLYMOL_HOME}/contrib/test/ring_replacement/6a.smi"), 0)
    self.assertEmpty(rr.process(m))

  def test_substructure_match(self):
    m = MolFromSmiles("O[1c]1cc[1c](OC)cc1")
    self.assertIsNotNone(m)

    rr = RingReplacement()
    self.assertTrue(rr.set_ring_atom_smarts("[1c]"))
    self.assertGreater(rr.read_replacement_rings("${LILLYMOL_HOME}/contrib/test/ring_replacement/6a.smi"), 0)
    products = rr.process(m)
    self.assertLen(products, 12)
    for p in products:
      self.assertEqual(p.number_isotopic_atoms(), 2)

  def test_support_requirement(self):
    m = MolFromSmiles("Oc1ccc(OC)cc1 start")
    self.assertIsNotNone(m)

    rr = RingReplacement()
    self.assertTrue(rr.set_ring_atom_smarts("c"))
    rr.set_min_support_requirement(100)
    self.assertGreater(rr.read_replacement_rings("${LILLYMOL_HOME}/contrib/test/ring_replacement/6a.smi"), 0)
    products = rr.process(m)
    self.assertLen(products, 7)
    for ndx, p in enumerate(products):
      self.assertEqual(p.number_isotopic_atoms(), 2, p.aromatic_smiles())

  def test_unique_only(self):
    m = MolFromSmiles("Oc1ccc(OC)cc1 start")
    self.assertIsNotNone(m)

    rr = RingReplacement()
    self.assertTrue(rr.set_ring_atom_smarts("c"))
    rr.set_unique_molecules_only(True)
    rr.set_min_support_requirement(100)
    self.assertGreater(rr.read_replacement_rings("${LILLYMOL_HOME}/contrib/test/ring_replacement/6a.smi"), 0)
    products1 = rr.process(m)
    self.assertLen(products1, 7)   # Seven generated products; the parent is not included
    products2 = rr.process(m)
    # All products were seen during the previous call.
    self.assertEmpty(products2)
    rr.clear_unique_molecule_cache()
    products3 = rr.process(m)
    self.assertLen(products3, 7)

if __name__ == '__main__':
  absltest.main()

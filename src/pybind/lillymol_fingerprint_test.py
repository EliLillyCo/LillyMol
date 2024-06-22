# Tests for fingerprinting.
import numpy as np

from absl import app
from absl import logging
from absl.testing import absltest


from lillymol import *
from lillymol_fingerprint import *

class TestECFingerprints(absltest.TestCase):
  def testMethane(self):
    mol = MolFromSmiles("C")
    self.assertEqual(mol.smiles(), "C")
    fp_creator = ECFingerprintCreator(1024);
    bits = fp_creator.fingerprint(mol)
    self.assertEqual(len(bits), 1024)
    self.assertEqual(bits.sum(), 1)

  def testOrderIndependent(self):
    mol = MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    fp_creator = ECFingerprintCreator(512);

    bits = fp_creator.fingerprint(mol)

    for i in range(100):
      smi = mol.random_smiles()
      m2 = MolFromSmiles(smi)
      self.assertEqual(mol.unique_smiles(), m2.unique_smiles());
      b = fp_creator.fingerprint(m2)

      self.assertEqual(bits.sum(), b.sum())
      self.assertTrue(np.array_equal(bits, b))

  def testIsSimilar(self):
    m1 = MolFromSmiles("CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4")
    m2 = MolFromSmiles("CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=NC=C4")

    fp_creator = ECFingerprintCreator(1024)
    bits1 = fp_creator.fingerprint(m1)
    bits2 = fp_creator.fingerprint(m2)
    self.assertFalse(np.array_equal(bits1, bits2))

    intersection = np.minimum(bits1, bits2)
    union = np.maximum(bits1, bits2)
    tanimoto = np.sum(intersection) / np.sum(union)
    self.assertAlmostEqual(tanimoto, 0.8531, places=3)

if __name__ == '__main__':
  absltest.main()

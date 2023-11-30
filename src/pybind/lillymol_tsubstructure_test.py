# Tests for tsubstructure object

from absl.testing import absltest

from lillymol import *
from lillymol_tsubstructure import *

class TestLillyMol(absltest.TestCase):
  def test_no_queries(self):
    ts = TSubstructure()
    self.assertEqual(ts.substructure_search("C methand"), 0)

  def test_substructure_search_single_query_single_molecule(self):
    ts = TSubstructure()
    self.assertTrue(ts.add_query_from_smarts("C carbon"))
    mol = MolFromSmiles("CC")
    self.assertTrue(ts.substructure_search(mol))

  def test_num_matches_single_query_single_molecule(self):
    ts = TSubstructure()
    self.assertTrue(ts.add_query_from_smarts("C carbon"))
    mol = MolFromSmiles("CC")
    self.assertEqual(ts.num_matches(mol), [2])

  def test_multiple_query_single_molecule(self):
    ts = TSubstructure()
    self.assertTrue(ts.add_query_from_smarts("C carbon"))
    self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
    mol = MolFromSmiles("CC")
    self.assertEqual(ts.num_matches(mol), [2, 0])

  def test_substructure_search_single_query_multiple_molecule(self):
    ts = TSubstructure()
    self.assertTrue(ts.add_query_from_smarts("C carbon"))
    smiles = ["C methane", "CC ethane", "CCC propane", "C1CC1 cyclopropane", "c1ccccc1 benzene"]
    mols = [MolFromSmiles(smi) for smi in smiles]
    self.assertEqual(ts.substructure_search(mols), [True, True, True, True, False])

  def test_num_matches_single_query_multiple_molecule(self):
    ts = TSubstructure()
    self.assertTrue(ts.add_query_from_smarts("C carbon"))
    smiles = ["C methane", "CC ethane", "CCC propane", "C1CC1 cyclopropane", "c1ccccc1 benzene"]
    mols = [MolFromSmiles(smi) for smi in smiles]
    self.assertEqual(ts.num_matches(mols), [[1], [2], [3], [3], [0]])

  def test_substructure_search_multiple_query_multiple_molecule(self):
    ts = TSubstructure()
    self.assertTrue(ts.add_query_from_smarts("C carbon"))
    self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
    smiles = ["C methane", "CC ethane", "N nitrogen", "O oxygem", "CN CN"]
    mols = [MolFromSmiles(smi) for smi in smiles]
    self.assertEqual(ts.substructure_search(mols), [True, True, True, False, True])

  def test_num_matches_multiple_query_multiple_molecule(self):
    ts = TSubstructure()
    self.assertTrue(ts.add_query_from_smarts("C carbon"))
    self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
    smiles = ["C methane", "CC ethane", "N nitrogen", "O oxygem", "CN CN"]
    mols = [MolFromSmiles(smi) for smi in smiles]
    self.assertEqual(ts.num_matches(mols), [[1, 0], [2, 0], [0, 1], [0, 0], [1, 1]])

  def test_must_match_all_queries(self):
    ts = TSubstructure()
    self.assertTrue(ts.add_query_from_smarts("C carbon"))
    self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
    mol = MolFromSmiles("C");
    self.assertTrue(ts.substructure_search(mol))
    ts.must_match_all_queries = True
    self.assertFalse(ts.substructure_search(mol))

  def test_unique_embeddings_only(self):
    ts = TSubstructure();
    mol = MolFromSmiles("CC(C)(C)c1c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c1C(C)(C)C")
    self.assertTrue(ts.add_query_from_smarts("CC(C)(C)c1c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c1C(C)(C)C"))
    self.assertEqual(ts.num_matches(mol), [559872])
    ts.set_unique_embeddings_only(True)
    self.assertEqual(ts.num_matches(mol), [1])

  def test_symmetry_equivalent(self):
    ts = TSubstructure();
    mol = MolFromSmiles("Cc1c(C)cccc1")
    self.assertTrue(ts.add_query_from_smarts("Cc"))
    self.assertEqual(ts.num_matches(mol), [2])
    ts.set_perceive_symmetry_equivalent_matches(False)
    self.assertEqual(ts.num_matches(mol), [1])

  def test_set_find_one_embedding_per_root_atom(self):
    ts = TSubstructure();
    mol = MolFromSmiles("CC(F)(F)F")
    self.assertTrue(ts.add_query_from_smarts("C(F)(F)F"))
    self.assertEqual(ts.num_matches(mol), [6])
    ts.set_find_one_embedding_per_root_atom(True)

  def test_set_max_matches_to_find(self):
    ts = TSubstructure();
    mol = MolFromSmiles("c1ccccc1")
    self.assertTrue(ts.add_query_from_smarts("c1ccccc1"))
    self.assertEqual(ts.num_matches(mol), [12])
    ts.set_max_matches_to_find(5)
    self.assertEqual(ts.num_matches(mol), [5])
    ts.set_max_matches_to_find(10)
    self.assertEqual(ts.num_matches(mol), [10])

  def test_reduce_to_largest_fragment(self):
    ts = TSubstructure();
    mol = MolFromSmiles("CC.C")
    self.assertTrue(ts.add_query_from_smarts("C"))
    self.assertEqual(ts.num_matches(mol), [3])
    ts.set_reduce_to_largest_fragment(True)
    self.assertEqual(ts.num_matches(mol), [2])

  def test_set_make_implicit_hydrogens_explicit(self):
    ts = TSubstructure();
    mol = MolFromSmiles("CC")
    self.assertTrue(ts.add_query_from_smarts("[#1]-C"))
    self.assertEqual(ts.num_matches(mol), [0])
    ts.set_make_implicit_hydrogens_explicit(True)
    self.assertEqual(ts.num_matches(mol), [6])

    ts.set_perceive_symmetry_equivalent_matches(False)
    self.assertEqual(ts.num_matches(mol), [1])

  def test_label_matched_atoms(self):
    ts = TSubstructure();
    mol = MolFromSmiles("Cc1ccccc1")
    self.assertTrue(ts.add_query_from_smarts("C"))
    self.assertTrue(ts.add_query_from_smarts("c"))
    ts.isotope = 4
    self.assertEqual(ts.label_matched_atoms(mol), 2)
    self.assertEqual(mol.unique_smiles(), "[4CH3][4c]1[4cH][4cH][4cH][4cH][4cH]1")

  def test_set_label_by_query_number(self):
    ts = TSubstructure();
    ts.set_label_by_query_number(True)
    self.assertTrue(ts.add_query_from_smarts("C"))
    self.assertTrue(ts.add_query_from_smarts("c"))
    self.assertTrue(ts.add_query_from_smarts("N"))
    self.assertTrue(ts.add_query_from_smarts("O"))
    self.assertTrue(ts.add_query_from_smarts("F"))

    mol = MolFromSmiles("Cc1cc(O)cc(F)c1")
    self.assertEqual(ts.label_matched_atoms(mol), 4)
    self.assertEqual(mol.unique_smiles(), "[5F][2c]1[2cH][2c]([1CH3])[2cH][2c]([4OH])[2cH]1")

if __name__ == '__main__':
  absltest.main()

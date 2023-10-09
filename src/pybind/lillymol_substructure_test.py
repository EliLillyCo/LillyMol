# Tests for substructure searching

import os
import tempfile

from absl import app
from absl import logging
from absl.testing import absltest

from google.protobuf import text_format

from Molecule_Lib import substructure_pb2

from lillymol import *
from lillymol_query import *

class TestLillyMolSubstructure(absltest.TestCase):
  def test_carbon(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("C methane"))
    qry = SubstructureQuery()
    self.assertTrue(qry.build_from_smarts("C"))
    self.assertEqual(qry.substructure_search(mol), 1)
    self.assertIn(qry, mol)

  def test_ethane_c(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("CC ethane"))
    qry = SubstructureQuery()
    self.assertTrue(qry.build_from_smarts("C"))
    self.assertEqual(qry.substructure_search(mol), 2)
    qry.set_perceive_symmetry_equivalent_matches(False)
    self.assertEqual(qry.substructure_search(mol), 1)

  def test_ethane_cc_embeddings_do_not_overlap(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("CC ethane"))
    qry = SubstructureQuery()
    self.assertTrue(qry.build_from_smarts("CC"))
    self.assertEqual(qry.substructure_search(mol), 2)
    qry.set_embeddings_do_not_overlap(True)
    self.assertEqual(qry.substructure_search(mol), 1)

  def test_ethane_cc_find_unique_embeddings_only(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("CC ethane"))
    qry = SubstructureQuery()
    self.assertTrue(qry.build_from_smarts("CC"))
    self.assertEqual(qry.substructure_search(mol), 2)
    qry.set_find_unique_embeddings_only(True)
    self.assertEqual(qry.substructure_search(mol), 1)

  def test_ethane_cc_find_one_embedding_per_atom(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("CC ethane"))
    qry = SubstructureQuery()
    self.assertTrue(qry.build_from_smarts("CC"))
    self.assertEqual(qry.substructure_search(mol), 2)
    qry.set_find_one_embedding_per_atom(True)
    self.assertEqual(qry.substructure_search(mol), 2)

  def test_results_returned(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("CC ethane"))
    qry = SubstructureQuery()
    self.assertTrue(qry.build_from_smarts("CC"))
    sresults = SubstructureResults()
    self.assertEqual(qry.substructure_search(mol, sresults), 2)

    self.assertEqual(sresults.number_embeddings(), 2)
    # Count the number of times each atoms is matched.
    matched = [0, 0]

    for embedding in sresults:
      self.assertLen(embedding, 2)
      for a in embedding:
        matched[a] += 1
    self.assertEqual(matched, [2, 2])

  def test_sresults_set_vector(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("Oc1c(O)c(O)c(O)c(O)c1O"))
    qry = SubstructureQuery()
    self.assertTrue(qry.build_from_smarts("Oc"))
    sresults = SubstructureResults()
    self.assertEqual(qry.substructure_search(mol, sresults), 6)
    v = sresults.each_embedding_set_vector(mol.natoms(), 2)
    self.assertEqual(v, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

  def test_matched_atoms_returned(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("Oc1ccccc1"))
    qry = SubstructureQuery()
    self.assertTrue(qry.build_from_smarts("Occ"))
    matches = qry.substructure_search_matches(mol)
    for match in matches:
      mol.set_isotopes(match, 1)
    self.assertEqual(mol.aromatic_smiles(), "[1OH][1c]1[1cH]ccc[1cH]1")

  def test_from_proto(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("Oc1ccccc1"))
    proto_string = """
query {
  smarts: "[OD1]-c:c"
  one_embedding_per_start_atom: true
}
"""
    proto = text_format.Parse(proto_string, substructure_pb2.SubstructureQuery())
    # This does not work yet. TODO:ianwatson
    # print(type(proto))
    # qry = SubstructureQuery()
    # self.assertTrue(qry.construct_from_proto(proto))

  def test_from_proto_file(self):
    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("Oc1ccccc1"))
    proto_string = """
query {
  smarts: "[OD1]-c:c"
  unique_embeddings_only: true
}
"""
    dir = absltest.TEST_TMPDIR.value
    tmpdir = tempfile.mkdtemp(dir=dir)
    fname = os.path.join(tmpdir, "qry.textproto")
    with open(fname, "w") as writer:
      writer.write(proto_string)
    qry = SubstructureQuery()
    self.assertTrue(qry.read_proto(fname))
    self.assertEqual(qry.substructure_search(mol), 2)
    self.assertIn(qry, mol)

if __name__ == '__main__':
  #app.run(absltest.main)
  absltest.main()

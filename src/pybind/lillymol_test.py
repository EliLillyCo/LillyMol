import copy
import math
import unittest

import numpy as np

from absl import app
from absl import logging
from absl.testing import absltest

from lillymol import *

CHIRAL_SMILES1 = [
"C(O)[C@@H](N)C CHEMBL1229871",
"O1[C@H](C1)CF CHEMBL501668",
"O1[C@H](C1)CBr CHEMBL504705",
"O1[C@H](C1)CCl CHEMBL448626",
"BrC[C@@H](C)O CHEMBL446288",
"C(Cl)[C@@H](C)Cl CHEMBL373466",
"SC[C@@H](N)C CHEMBL37279",
"OC[C@H](N)CC CHEMBL3184640",
"C1(=O)[C@H](N)CO1 CHEMBL2219717",
"ClC[C@H](O)CO CHEMBL1794186",
]
CHIRAL_SMILES2 = [
"O[C@H]1[C@H](N)CCC1 CHEMBL2374489",
"O[C@H]1[C@H](N)CCC1 CHEMBL2375114",
"N1C[C@H](O)[C@@H](O)C1 CHEMBL2335511",
"C1[C@@H](O)[C@@H](O)CN1 CHEMBL2207396",
"C1[C@@H](N)[C@H]1C(=O)O CHEMBL403157",
"C1NC[C@@H](O)[C@@H]1O CHEMBL396701",
"C1C[C@@H](O)[C@@H](O)C1 CHEMBL399324",
"C1=NC[C@@H](O)[C@H]1O CHEMBL389969",
"O1[C@H](C)[C@@H]1C(=O)O CHEMBL370643",
"O1C[C@@H](O)[C@H](O)C1 CHEMBL350524",
]

class TestLillyMol(absltest.TestCase):
  def test_empty_molecule(self):
    m = Molecule()
    self.assertTrue(m.empty())
    self.assertEmpty(m)
    self.assertEqual(m.natoms(), 0)
    self.assertEqual(m.name(), "")

  def test_methane(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C methane"))
    self.assertEqual(m.name(), "methane")
    self.assertEqual(m.natoms(), 1)
    self.assertEqual(m.natoms(6), 1)
    self.assertEqual(m.natoms(7), 0)
    self.assertEqual(m.nedges(), 0)
    self.assertEqual(m.nrings(), 0)
    self.assertFalse(m.is_ring_atom(0))
    self.assertEqual(m.number_fragments(), 1)
    self.assertEqual(m.atomic_symbol(0), "C")
    self.assertEqual(m.atomic_number(0), 6)
    self.assertEqual(m.molecular_formula(), "CH4")
    self.assertEqual(m.smiles(), "C")
    self.assertEqual(m.unique_smiles(), "C")
    self.assertEmpty(m.rings())
    self.assertEqual(m.hcount(0), 4)
    self.assertEqual(m.implicit_hydrogens(0), 4)
    self.assertEqual(m.explicit_hydrogens(0), 0)
    self.assertEqual(m.isotope(0), 0)
    self.assertTrue(m.valence_ok())
    self.assertEqual(m.highest_coordinate_dimensionality(), 0)

  def test_copy_constructor(self):
    m1 = Molecule()
    self.assertTrue(m1.build_from_smiles("C"))
    m2 = Molecule(m1)
    self.assertEqual(m2.smiles(), "C")
    m1.add_atom(6)
    self.assertEqual(m1.smiles(), "C.C")
    self.assertEqual(m2.smiles(), "C")

  def test_copy_constructor_with_name(self):
    m = Molecule()
    set_copy_name_in_molecule_copy_constructor(True)
    self.assertTrue(m.build_from_smiles("C methane"))
    m2 = Molecule(m)
    set_copy_name_in_molecule_copy_constructor(False)
    self.assertEqual(m2.name(), "methane")

  def test_copy_operation(self):
    m1 = Molecule()
    self.assertTrue(m1.build_from_smiles("C"))
    m2 = copy.copy(m1)
    self.assertEqual(m2.smiles(), "C")
    m2.add_atom(6)
    self.assertEqual(m1.smiles(), "C")
    self.assertEqual(m2.smiles(), "C.C")

  def test_hydrogen_related(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C([H])([H])([H])N([H])C([H])([H])C([H])([H])C([H])(C1=C([H])C(=C([H])C(=C1[H])[H])[H])OC1=C([H])C(=C(C(=C1[H])[H])C(F)(F)F)[H] prozac"))
    self.assertEqual(m.natoms(), 22 + 18)
    self.assertEqual(m.natoms(1), 18)

    self.assertEqual(m.atomic_number(0), 6)
    self.assertEqual(m.ncon(0), 4)
    self.assertEqual(m.ncon(1), 1)
    self.assertAlmostEqual(m.amw(), 309.33, 2)
    self.assertAlmostEqual(m.exact_mass(), 309.13404868, 6)
    m.remove_all(1)
    self.assertAlmostEqual(m.amw(), 309.33, 2)
    self.assertAlmostEqual(m.exact_mass(), 309.13404868, 6)
    m.make_implicit_hydrogens_explicit()
    self.assertEqual(m.natoms(), 22 + 18)
    self.assertEqual(m.natoms(1), 18)
    self.assertAlmostEqual(m.amw(), 309.33, 2)
    self.assertAlmostEqual(m.exact_mass(), 309.13404868, 6)

  def test_isotopes(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("[1CH4]"))
    self.assertTrue(m.valence_ok())
    self.assertEqual(m.isotope(0), 1)
    self.assertEqual(m.number_isotopic_atoms(), 1)
    m.set_isotope(0, 2)
    self.assertEqual(m.smiles(), "[2CH4]")
    m.set_isotope(0, 108)
    self.assertEqual(m.smiles(), "[108CH4]")
    m.remove_isotopes()
    self.assertEqual(m.smiles(), "C")

    self.assertTrue(m.build_from_smiles("C[1C][2C][3C]C"))
    self.assertFalse(m.valence_ok())
    self.assertEqual(m.natoms(1), 0)
    self.assertEqual(m.hcount(0), 3)
    self.assertEqual(m.hcount(1), 0)
    self.assertEqual(m.hcount(2), 0)
    self.assertEqual(m.hcount(3), 0)
    self.assertEqual(m.hcount(4), 3)
    for i in range(m.natoms()):
      m.unset_all_implicit_hydrogen_information(i)
    self.assertTrue(m.valence_ok())
    self.assertEqual(m.hcount(1), 2)
    self.assertEqual(m.hcount(2), 2)
    self.assertEqual(m.hcount(3), 2)

    s = Set_of_Atoms([1, 2, 3])
    m.set_isotopes(s, 5)
    self.assertEqual(m.isotope(1), 5)
    self.assertEqual(m.isotope(2), 5)
    self.assertEqual(m.isotope(3), 5)

    m.remove_isotopes()
    self.assertEqual(m.smiles(), "CCCCC")

    self.assertTrue(m.build_from_smiles("C[1C][2C][3C]C"))
    self.assertFalse(m.valence_ok())
    m.remove_hydrogens_known_flag_to_fix_valence_errors()
    self.assertTrue(m.valence_ok())

    m = Molecule()
    self.assertTrue(m.build_from_smiles("C[1C][2C][3C][1C]"))
    self.assertEqual(m.first_atom_with_isotope(2), 2)
    self.assertEqual(m.first_atom_with_isotope(1), 1)
    self.assertEqual(m.first_atom_with_isotope(7), -1)


  def test_fragment_related(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C1CC1C1CC1"))
    self.assertEqual(m.natoms(), 6)
    self.assertEqual(m.nrings(), 2)
    self.assertEqual(m.number_fragments(), 1)
    self.assertTrue(m.in_same_ring(0, 1))
    self.assertTrue(m.in_same_ring(1, 2))
    self.assertFalse(m.in_same_ring(2, 3))
    self.assertTrue(m.in_same_ring(3, 4))
    self.assertTrue(m.in_same_ring(4, 5))
    self.assertFalse(m.in_same_ring(2, 3))

    m.remove_bond_between_atoms(2, 3)
    self.assertEqual(m.natoms(), 6)
    self.assertEqual(m.nrings(), 2)
    self.assertEqual(m.number_fragments(), 2)
    self.assertTrue(m.in_same_ring(0, 1))
    self.assertTrue(m.in_same_ring(1, 2))
    self.assertFalse(m.in_same_ring(2, 3))
    self.assertTrue(m.in_same_ring(3, 4))
    self.assertTrue(m.in_same_ring(4, 5))
    self.assertNotEqual(m.fused_system_identifier(2), m.fused_system_identifier(3))
    self.assertEqual(m.atoms_in_fragment(0), 3)
    self.assertEqual(m.atoms_in_fragment(1), 3)
    self.assertEqual(m.fragment_membership(0), m.fragment_membership(1))
    self.assertEqual(m.fragment_membership(1), m.fragment_membership(2))
    self.assertEqual(m.fragment_membership(3), m.fragment_membership(4))
    self.assertEqual(m.fragment_membership(4), m.fragment_membership(5))
    # The actual numbers should not be depended on, this test is likely
    # unstable
    self.assertEqual(m.label_atoms_by_ring_system(), [1, 1, 1, 2, 2, 2])

    # Restore to where we started
    m.add_bond(2, 3, BondType.SINGLE_BOND)
    self.assertEqual(m.number_fragments(), 1)
    self.assertEqual(m.smiles(), "C1CC1C1CC1")
    self.assertEqual(m.number_ring_systems(), 2)
    self.assertEqual(m.label_atoms_by_ring_system(), [1, 1, 1, 2, 2, 2])

  def test_atoms_in_fragment(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CCCC.c1ccccc1"))
    self.assertEqual(m.atoms_in_fragment(0), 4);
    self.assertEqual(m.atoms_in_fragment(1), 6);

  def test_aspirin(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O aspirin"))
    self.assertEqual(m.molecular_formula(), "C9H8O4")
    self.assertAlmostEqual(m.amw(), 180.16, 2)
    self.assertAlmostEqual(m.exact_mass(), 180.04225873, 6)
    self.assertEqual(m.natoms(), 13)
    self.assertEqual(m.nedges(), 13)
    self.assertEqual(m.nrings(), 1)
    self.assertEqual(m.aromatic_ring_count(), 1)
    self.assertEqual(m.aromatic_atom_count(), 6)
    self.assertEqual(m.connections(0), [1])
    self.assertEqual(m.number_fragments(), 1)
    self.assertEqual(m.connections(1), [0, 2, 3])
    self.assertTrue(m.organic_only())
    self.assertTrue(m.is_aromatic(4))

  def test_largest_fragment(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C.CC.CCC"))
    self.assertEqual(m.natoms(), 6)
    self.assertEqual(m.number_fragments(), 3)
    m.reduce_to_largest_fragment()
    self.assertEqual(m.smiles(), "CCC")

  def test_largest_fragment_carefully(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("O=N(=O)C1=CC(=C(O)C(=C1)N(=O)=O)N(=O)=O.CCSC(=N)N CHEMBL4531203"))
    m.reduce_to_largest_fragment_carefully()
    self.assertEqual(m.smiles(), "O=N(=O)C1=CC(=C(O)C(=C1)N(=O)=O)N(=O)=O")

    self.assertTrue(m.build_from_smiles("O(C(=O)C(C1=NC(C)(C)CC2=CC=CC=C12)CC)CC.OC1=C(N(=O)=O)C=C(N(=O)=O)C=C1N(=O)=O CHEMBL1352385"))
    m.reduce_to_largest_fragment_carefully()
    self.assertEqual(m.smiles(), "O(C(=O)C(C1=NC(C)(C)CC2=CC=CC=C12)CC)CC")

    self.assertTrue(m.build_from_smiles("N1=C(NC(C)C1)CC1=CC=CC=C1.C1(=CC(=CC(=C1O)N(=O)=O)N(=O)=O)N(=O)=O CHEMBL609421"))
    m.reduce_to_largest_fragment_carefully()
    self.assertEqual(m.smiles(), "N1=C(NC(C)C1)CC1=CC=CC=C1")

  def test_fragment_membership(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C.CC.CCC.CCCC"))
    self.assertEqual(m.get_fragment_membership(), [0, 1, 1, 2, 2, 2, 3, 3, 3, 3])

  def test_delete_fragment(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC.CCC.CC"))
    # Fragment numbering is not guaranteed. Better to programatically determine
    # an atom in the fragment and remove the fragment containing that atom.
    m.delete_fragment(1)
    self.assertEqual(m.smiles(), "CC.CC")

    # Remove the fragment containing a nitrogen atom
    self.assertTrue(m.build_from_smiles("CC.CNC.CC"))
    to_remove = -1
    for ndx, atom in enumerate(m):
      if atom.atomic_number() == 7:
        to_remove = m.fragment_membership(ndx)
        break
    self.assertGreaterEqual(to_remove, 0)
    m.delete_fragment(to_remove)
    self.assertEqual(m.smiles(), "CC.CC")

  def test_saturated(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC=CC#CCc1ccccc1"))
    self.assertTrue(m.saturated(0))
    self.assertFalse(m.saturated(1))
    self.assertFalse(m.saturated(2))
    self.assertFalse(m.saturated(3))
    self.assertFalse(m.saturated(4))
    self.assertFalse(m.saturated(6))

  def test_unsaturation(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC=CC#CCc1ccccc1"))
    self.assertEqual(m.unsaturation(0), 0)
    self.assertEqual(m.unsaturation(1), 1)
    self.assertEqual(m.unsaturation(2), 1)
    self.assertEqual(m.unsaturation(3), 2)
    self.assertEqual(m.unsaturation(4), 2)
    self.assertEqual(m.unsaturation(5), 0)
    self.assertEqual(m.unsaturation(6), 1)

  def test_remove_atom(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CCC"))
    m.remove_atom(1)
    self.assertEqual(m.smiles(), "C.C")

    m.add_atom(7)
    m.add_bond(0, 2, BondType.SINGLE_BOND)
    m.add_bond(1, 2, BondType.SINGLE_BOND)
    self.assertEqual(m.smiles(), "CNC")

  def test_set_bond_type_between_atoms(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC"))
    m.set_bond_type_between_atoms(0, 1, BondType.DOUBLE_BOND)
    self.assertEqual(m.smiles(), "C=C")
    m.set_bond_type_between_atoms(0, 1, BondType.TRIPLE_BOND)
    self.assertEqual(m.smiles(), "C#C")

  def test_set_atomic_number(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC"))
    m.set_atomic_number(1, 103)
    self.assertEqual(m.smiles(), "C[Lr]")

  def test_set_auto_create_new_elements(self):
    m = Molecule()
    set_auto_create_new_elements(True)
    self.assertTrue(m.build_from_smiles("[Th][Eq][U]IC[K]BrO[W]NFO[Xj][Um]PSO[Ve][Rt][He][La][Zy][D]O[G]"))
    self.assertEqual(m.natoms(), 25)
    self.assertTrue("Br" in m)
    self.assertTrue("Um" in m)
    set_auto_create_new_elements(False)

  def test_set_atomic_symbols_can_have_arbitrary_length(self):
    m = Molecule()
    set_atomic_symbols_can_have_arbitrary_length(True)
    set_auto_create_new_elements(True)
    self.assertTrue(m.build_from_smiles("[Ala][Glu]"))
    self.assertEqual(m.natoms(), 2)
    self.assertEqual(m.atomic_symbol(0), "Ala")
    self.assertEqual(m.atomic_symbol(1), "Glu")
    self.assertTrue("Ala" in m)
    set_auto_create_new_elements(False)
    set_atomic_symbols_can_have_arbitrary_length(False)

  def test_fused_to(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C1C2CC12"))
    self.assertEqual(m.nrings(), 2)
    self.assertLen(m.ring(0), 3)
    self.assertLen(m.ring(1), 3)
    r0 = m.ring(0)
    r1 = m.ring(1)
    self.assertEqual(r0.fused_system_identifier(), r1.fused_system_identifier())
    self.assertTrue(r0.is_fused())
    self.assertTrue(r1.is_fused())
    self.assertTrue(r0.is_fused_to(r1))
    self.assertTrue(r1.is_fused_to(r0))

    self.assertEqual(r0.fragment_membership(), r1.fragment_membership())

    self.assertEqual(r0.largest_number_of_bonds_shared_with_another_ring(), 1)
    self.assertEqual(r1.largest_number_of_bonds_shared_with_another_ring(), 1)

  def test_add(self):
    m1 = Molecule()
    m2 = Molecule()
    self.assertTrue(m1.build_from_smiles("C"))
    self.assertTrue(m2.build_from_smiles("C"))
    m1.add(m2)
    self.assertEqual(m1.smiles(), "C.C")
    # The donor is unchanged
    self.assertEqual(m2.smiles(), "C")

  def test_cubane(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C12C3C4C1C5C2C3C45"))
    self.assertEqual(m.nrings(), 5)
    self.assertEqual(m.non_sssr_rings(), 1)
    for i in range(m.natoms()):
      self.assertEqual(m.ring_bond_count(i), 3)

    for ring in m.rings():
      self.assertLen(ring, 4)

    # The order of ring finding is not guaranteed
    self.assertEqual(m.get_ring_membership(), [3, 3, 3, 3, 2, 2, 2, 2])

    self.assertEqual(m.fused_system_identifier(0), m.fused_system_identifier(1))

  def test_distance_matrix(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C1CCC1"))
    self.assertEqual(m.bonds_between(0, 1), 1)
    self.assertEqual(m.bonds_between(0, 2), 2)

    self.assertTrue(m.build_from_smiles("N(CC1=CC=C(OCCCC2=CC=CC=C2)C=C1)(CC1=CC=C(OCCCC2=CC=CC=C2)C=C1)CCCCN CHEMBL349114"))
    self.assertEqual(m.fragment_membership(30), m.fragment_membership(39))
    self.assertEqual(m.bonds_between(30, 39), 18)
    self.assertEqual(m.longest_path(), 26)
    self.assertEqual(m.most_distant_pair(), (13, 30))
    self.assertEqual(m.bonds_between(13, 30), 26)

  def test_atom_numbers(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C[CH2:1][CH2:2]C"))
    self.assertTrue(m.valence_ok())
    self.assertEqual(m.atom_map_number(0), 0)
    self.assertEqual(m.atom_map_number(1), 1)
    self.assertEqual(m.atom_map_number(2), 2)
    self.assertEqual(m.atom_map_number(3), 0)

    m.set_atom_map_number(3, 3)
    self.assertEqual(m.smiles(), "C[CH2:1][CH2:2][CH3:3]")

    self.assertEqual(m.atom_with_atom_map_number(3), 3)

    m.reset_atom_map_numbers()
    self.assertEqual(m.smiles(), "CCCC")

  def test_sort_atoms(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CO.N"))
    self.assertEqual(m.atomic_number(0), 6)
    self.assertEqual(m.atomic_number(1), 8)
    self.assertEqual(m.atomic_number(2), 7)
    order = [1, 2, 0]
    m.sort_atoms(order)
    self.assertEqual(m.atomic_number(0), 8)
    self.assertEqual(m.atomic_number(1), 6)
    self.assertEqual(m.atomic_number(2), 7)
    self.assertEqual(m.smiles(), "OC.N")


  def test_create_components(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C.CC.CCC.CCCC"))
    components = m.create_components()
    self.assertLen(components, 4)
    # The order of the component molecules is not guaranteed.
    self.assertEqual(components[0].natoms(), 1)
    self.assertEqual(components[1].natoms(), 2)
    self.assertEqual(components[2].natoms(), 3)
    self.assertEqual(components[3].natoms(), 4)

  def test_random_smiles(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O aspirin"))
    initial_smiles = m.unique_smiles()
    seen = set()
    for i in range(10):
      smiles = m.random_smiles()
      if smiles in seen:
        continue
      seen.add(smiles)
      m2 = Molecule()
      self.assertTrue(m2.build_from_smiles(smiles))
      self.assertEqual(m2.unique_smiles(), initial_smiles)

  def test_smiles_starting_atom(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("OCN"))
    self.assertEqual(m.smiles(), "OCN")
    self.assertEqual(m.smiles_starting_with_atom(1), "C(O)N")

  def test_isotopically_labelled_smiles(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("OCN"))
    self.assertEqual(m.isotopically_labelled_smiles(), "O[1CH2][2NH2]")

  def test_symmetry(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("NC1=CC=C(C=C1)C(F)(F)F CHEMBL1162294"))
    self.assertEqual(m.number_symmetry_classes(), 7)
    self.assertEqual(m.symmetry_class(2), m.symmetry_class(6))
    self.assertEqual(m.symmetry_class(3), m.symmetry_class(5))

    self.assertEqual(m.symmetry_class(8), m.symmetry_class(9))
    self.assertEqual(m.symmetry_class(9), m.symmetry_class(10))
    self.assertEqual(m.symmetry_equivalents(8), [9, 10])
    self.assertEmpty(m.symmetry_equivalents(0))
    self.assertEmpty(m.symmetry_equivalents(4))
    self.assertEmpty(m.symmetry_equivalents(7))

  def test_organic(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("OCN"))
    self.assertTrue(m.organic_only())
    self.assertTrue(m.build_from_smiles("OCNB"))
    self.assertFalse(m.organic_only())

  def test_chop(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("OCN"))
    m.chop(1)
    self.assertEqual(m.smiles(), "OC")

  def test_remove_bond(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("OCN"))
    m.remove_bond_between_atoms(0, 1)
    self.assertEqual(m.smiles(), "O.CN")

  def test_remove_bonds_to_atom(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("OCN"))
    m.remove_bonds_to_atom(1);
    self.assertEqual(m.smiles(), "O.C.N")

  def test_remove_all_atomic_number(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("OCN"))
    m.remove_all(8)
    self.assertEqual(m.smiles(), "CN")


  #      C
  #    / | \
  # C-C  |  C-C
  #    \ | /
  #      C
  #
  def test_ring_related(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC1C2C(C)C12"))
    self.assertEqual(m.natoms(), 6)
    self.assertEqual(m.natoms(6), 6)

    self.assertEqual(m.ncon(0), 1)
    self.assertEqual(m.ncon(1), 3)
    self.assertEqual(m.ncon(2), 3)
    self.assertEqual(m.ncon(3), 3)
    self.assertEqual(m.ncon(4), 1)
    self.assertEqual(m.ncon(5), 3)

    self.assertEqual(m.nrings(), 2)
    self.assertEqual(m.non_sssr_rings(), 0)
    self.assertEqual(m.nrings(0), 0)
    self.assertEqual(m.nrings(1), 1)
    self.assertEqual(m.nrings(2), 2)
    self.assertEqual(m.nrings(3), 1)
    self.assertEqual(m.nrings(4), 0)
    self.assertEqual(m.nrings(5), 2)

    self.assertEqual(m.ring_bond_count(0), 0)
    self.assertEqual(m.ring_bond_count(1), 2)
    self.assertEqual(m.ring_bond_count(2), 3)
    self.assertEqual(m.ring_bond_count(3), 2)
    self.assertEqual(m.ring_bond_count(4), 0)
    self.assertEqual(m.ring_bond_count(5), 3)

    for i in range(m.natoms()):
      self.assertEqual(m.attached_heteroatom_count(i), 0)

    self.assertFalse(m.is_ring_atom(0))
    self.assertTrue(m.is_ring_atom(1))

    self.assertFalse(m.in_ring_of_given_size(0, 3))
    self.assertTrue(m.in_ring_of_given_size(1, 3))
    self.assertFalse(m.in_ring_of_given_size(1, 4))
    self.assertTrue(m.in_ring_of_given_size(2, 3))
    self.assertTrue(m.in_ring_of_given_size(3, 3))
    self.assertTrue(m.in_ring_of_given_size(5, 3))

    self.assertTrue(m.in_same_ring(1, 2))
    self.assertTrue(m.in_same_ring(1, 5))
    self.assertTrue(m.in_same_ring(3, 2))
    self.assertTrue(m.in_same_ring(3, 5))
    self.assertFalse(m.in_same_ring(1, 3))

    self.assertEqual(m.fused_system_identifier(1), m.fused_system_identifier(2))
    self.assertEqual(m.fused_system_identifier(1), m.fused_system_identifier(3))
    self.assertEqual(m.fused_system_identifier(1), m.fused_system_identifier(5))
    self.assertEqual(m.fused_system_identifier(1), m.fused_system_identifier(5))

    self.assertEqual(m.fused_system_size(1), 2)
    self.assertEqual(m.fused_system_size(0), 0)

    self.assertLen(m.rings(), 2)
    # The order of the rings should not be relied upon, so this test is
    # possibly fragile.
    self.assertCountEqual(m.ring(0), [1, 2, 5])
    self.assertCountEqual(m.ring(1), [2, 3, 5])
    for r in m.rings():
      self.assertLen(r, 3)

    self.assertEqual(m.largest_ring_size(), 3)
    self.assertEqual(m.number_ring_systems(), 1)
    self.assertFalse(m.is_spiro_fused(1))

  def test_atom_iterator(self):
    m = Molecule()
    #self.assertTrue(m.build_from_smiles("N1C(=O)CC1SC1=CC=C(F)C(=C1)Cl CHEMBL3394616")
    self.assertTrue(m.build_from_smiles("C1CC1"))
    self.assertTrue(m.valence_ok())
    # Iterate through all atoms in the molecule.
    for ndx, atom in enumerate(m):
      self.assertEqual(atom.atomic_number(), 6)
      self.assertEqual(m.ring_bond_count(ndx), 2)

    self.assertTrue(m.build_from_smiles("CC(N)(O)F"))
    # Loop through bonds connected to atom 1. Check with `m` that
    # the atom is bonded to 1.
    for bond in m[1]: 
      atom = bond.other(1)
      self.assertTrue(m.are_bonded(atom, 1))

  def test_iterate_bonds(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("O=C1N(C(=O)N2[C@@]1([C@@H]1[C@H]([C@H]2C2=CC=CC(=C2)C)C(=O)N(C1=O)CC1=CC=CC=C1)CC)C1=CC=CC=C1 CHEMBL1434237"))
    # Count number of aromatic bonds
    m.compute_aromaticity_if_needed()
    aromatic_bonds = 0
    for bond in m.bonds():
      if bond.is_aromatic():
        aromatic_bonds += 1
    self.assertEqual(aromatic_bonds, 18)

  def test_formal_charge(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C[N+H3] methylamine"))
    self.assertEqual(m.number_formal_charges(), 1)
    self.assertEqual(m.formal_charge(0), 0)
    self.assertEqual(m.formal_charge(1), 1)
    self.assertEqual(m.net_formal_charge(), 1)

    m.set_formal_charge(1, 0)
    self.assertEqual(m.smiles(), "CN")

  def build_benzene(self):
    m = Molecule()
    for i in range(6):
      m.add(6)

    self.assertEqual(m.natoms(), 6)
    self.assertEqual(m.natoms(6), 6)

    m.add_bond(0, 1, BondType.SINGLE_BOND)
    m.add_bond(1, 2, BondType.DOUBLE_BOND)
    m.add_bond(2, 3, BondType.SINGLE_BOND)
    m.add_bond(3, 4, BondType.DOUBLE_BOND)
    m.add_bond(4, 5, BondType.SINGLE_BOND)
    m.add_bond(5, 0, BondType.DOUBLE_BOND)
    self.assertEqual(m.aromatic_ring_count(), 1)

  def test_atom_iterator_and_valence(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("[1CH3]C(=[2CH2])#[3CH]"))
    self.assertFalse(m.valence_ok())

    for bond in m[1]:
      other = bond.other(1)
      if bond.is_single_bond():
        self.assertEqual(m.isotope(other), 1)
      elif bond.is_double_bond():
        self.assertEqual(m.isotope(other), 2)
      elif bond.is_triple_bond():
        self.assertEqual(m.isotope(other), 3)

  def test_atom_subtraction(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC{{1,1,1}}"))
    self.assertAlmostEqual(m[0] - m[1], math.sqrt(3.0))

  # Do a 'substructure search'. Generally do not do this, do an
  # actual substructure search, but this might be illustrative.
  # Identify a O=a atom pair.
  # There are many ways this could be done, and which will be most
  # efficient will depend on the molecules being examined.
  def test_find_exocyclic_bond(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("N1NC(=O)C=C1 CHEMBL4227850"))
    self.assertEqual(m.aromatic_ring_count(), 1)

    # Start search with singly bonded O atoms
    found_match = 0
    for ndx, atom in enumerate(m):
      if atom.ncon() != 1:
        continue
      if atom.atomic_number() != 8:
        continue
      # Bond to the first neighbour
      bond = atom[0]
      if not bond.is_double_bond():
        continue
      other = bond.other(ndx)
      if m.is_aromatic(other):
        found_match += 1

    self.assertGreater(found_match, 0)

    # Start search by looking at atoms in aromatic rings
    found_match = 0
    m.compute_aromaticity_if_needed()
    for ring in m.rings():
      if not ring.is_aromatic():
        continue
      for atom_number in ring:
        atom = m[atom_number]
        if atom.ncon() == 2:
          continue
        for bond in atom:
          if not bond.is_double_bond():
            continue
          other = bond.other(atom_number)
          if m.ncon(other) != 1:
            continue
          if m.atomic_number(other) == 8:
            found_match += 1
            break
    self.assertGreater(found_match, 0)

    # Start search with 3 connected aromatic atoms
    found_match = 0
    for ndx, atom in enumerate(m):
      if atom.ncon() != 3:
        continue
      if not m.is_aromatic(ndx):
        continue
      for bond in atom:
        if not bond.is_double_bond():
          continue
        other = bond.other(ndx)
        if m.ncon(other) != 1:
          continue
        if m.atomic_number(other) != 8:
          continue
        found_match += 1
    self.assertGreater(found_match, 0)

    # On 50k random Chembl molecules the times for the methods are (seconds)
    # 1 5.18
    # 2 9.14
    # 3 8.47
    # In this case, the singly bonded oxygen is the best place to start
    # since they are comparatively rare.

  def test_scaffold(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C"))
    m.to_scaffold()
    self.assertEqual(m.smiles(), "C")

    self.assertTrue(m.build_from_smiles("O1N=C(C(=O)N2CCCC2)C=C1COC1=CC=C2N=CC=CC2=C1 CHEMBL1589003"))
    m.to_scaffold()
    self.assertEqual(m.smiles(), "O1N=C(CN2CCCC2)C=C1COC1=CC=C2N=CC=CC2=C1")

    self.assertTrue(m.build_from_smiles("O=C(N(C1=CC=C(C)C=C1)CC(=O)NCCOC)CCC(=O)NC1=CC=CC=N1 CHEMBL1576099"))
    m.to_scaffold()
    self.assertEqual(m.smiles(), "C(NC1=CC=CC=C1)CCCNC1=CC=CC=N1")

    self.assertTrue(m.build_from_smiles("O=C1N(C(=O)C2=C1C(=CC=C2)N(=O)=O)CC(=O)N1CC2=CC=CC=C2CC1 CHEMBL2134451"))
    m.to_scaffold()
    self.assertEqual(m.smiles(), "C1N(CC2=C1C=CC=C2)CCN1CC2=CC=CC=C2CC1")

    self.assertTrue(m.build_from_smiles("O=C(C1=CC=CN1CC(=O)NCC1N(CCC1)CC)C1=CC=CC=C1C CHEMBL1404612"))
    m.to_scaffold()
    self.assertEqual(m.smiles(), "C(C1=CC=CN1CCNCC1NCCC1)C1=CC=CC=C1")

    self.assertTrue(m.build_from_smiles("O=C(N1[C@H](C(=O)NC2C3=CC=CC=C3CCC2)CCC1)[C@@H](NC(=O)[C@H](C)NC)CC(=O)O CHEMBL1570483"))
    m.to_scaffold()
    self.assertEqual(m.smiles(), "N1[C@H](CNC2C3=CC=CC=C3CCC2)CCC1")

  def test_coords(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C{{0,0,0}}C{{1,1,1}}C{{2,0,0}}"))
    self.assertEqual(m.highest_coordinate_dimensionality(), 3)
    self.assertAlmostEqual(m.x(0), 0.0)
    self.assertAlmostEqual(m.y(0), 0.0)
    self.assertAlmostEqual(m.z(0), 0.0)

    self.assertAlmostEqual(m.x(1), 1.0)
    self.assertAlmostEqual(m.y(1), 1.0)
    self.assertAlmostEqual(m.z(1), 1.0)

    self.assertAlmostEqual(m.distance_between_atoms(0, 1), math.sqrt(3.0))
    self.assertAlmostEqual(m.distance_between_atoms(0, 2), 2.0)

    self.assertAlmostEqual(m.bond_angle(0, 1, 2), 1.230959, 4)

  def test_getxyz(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C{{0,0,0}}C{{1,1,1}}C{{2,0,0}}"))
    coords = m.get_coordinates()
    self.assertEqual(len(coords), m.natoms() * 3)
    self.assertAlmostEqual(coords[0], 0.0)
    self.assertAlmostEqual(coords[1], 0.0)
    self.assertAlmostEqual(coords[2], 0.0)

    self.assertAlmostEqual(coords[3], 1.0)
    self.assertAlmostEqual(coords[4], 1.0)
    self.assertAlmostEqual(coords[5], 1.0)

    self.assertAlmostEqual(coords[6], 2.0)
    self.assertAlmostEqual(coords[7], 0.0)
    self.assertAlmostEqual(coords[8], 0.0)

  def test_setxyz(self):
    coords = np.arange(0, 6, dtype=float)
    self.assertEqual(len(coords), 6)

    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC"))
    self.assertAlmostEqual(m.x(0), 0.0)
    self.assertAlmostEqual(m.y(0), 0.0)
    self.assertAlmostEqual(m.z(0), 0.0)
    self.assertAlmostEqual(m.x(1), 0.0)
    self.assertAlmostEqual(m.y(1), 0.0)
    self.assertAlmostEqual(m.z(1), 0.0)

    m.set_coordinates(coords)
    for i in range(0, 2):
      self.assertAlmostEqual(m.x(i), i * 3 + 0)
      self.assertAlmostEqual(m.y(i), i * 3 + 1)
      self.assertAlmostEqual(m.z(i), i * 3 + 2)

  def test_dihedral_scan(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C{{-2,1,0}}C{{-1,0,0}}C{{0,0,0}}C{{1,1,0}}"))
    bump_check = 0.0 
    angle = 45.0
    coords = m.dihedral_scan(1, 2, angle, bump_check)
    self.assertEqual(len(coords), 7)

    expected = [-45.0, -90.0, -135.0, 180, 135, 90, 45]
    for ndx, c in enumerate(coords):
      m.set_coordinates(c)
      found = m.signed_dihedral_angle(0, 1, 2, 3)
      self.assertAlmostEqual(found * 180.0 / 3.14159265, expected[ndx], places=4)
      
  def test_rule_of_five(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("COC1=CC=C(C=C1)N2C3=C(CCN(C3=O)C4=CC=C(C=C4)N5CCCCC5=O)C(=N2)C(=O)N Eliquis"))
    donor = 0
    acceptor = 0

    # no consideration of charged atoms.
    for ndx,atom in enumerate(m):
      # Intercept the most common case first
      if atom.atomic_number() == 6:
        continue

      hcount = m.hcount(ndx)
      if atom.atomic_number() == 7:
        if hcount == 0:
          acceptor += 1
        else:
          donor += hcount
      elif atom.atomic_number() == 8:
        if hcount == 0:
          acceptor += 1
        else:
          donor += 1

    self.assertEqual(donor, 2)
    self.assertEqual(acceptor, 8)


  def test_number_chiral_smiles1(self):
    mols = [LillyMolFromSmiles(s) for s in CHIRAL_SMILES1];
    for i,m in enumerate(mols):
      self.assertEqual(m.number_chiral_centres(), 1)

  def test_number_chiral_smiles2(self):
    mols = [LillyMolFromSmiles(s) for s in CHIRAL_SMILES2];
    for i,m in enumerate(mols):
      self.assertEqual(m.number_chiral_centres(), 2)

  def test_chiral_implicit_hydrogen(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C[C@H](N)F"))
    self.assertEqual(m.number_chiral_centres(), 1)
    self.assertIsNone(m.chiral_centre_at_atom(0))
    c = m.chiral_centre_at_atom(1)
    self.assertIsNotNone(c)
    # Which atom gets assigned to which position is unpredictable.
    self.assertEqual(c.top_front(), 0)
    self.assertEqual(c.left_down(), 2)
    self.assertEqual(c.right_down(), 3)
    # TODO:ianwatson implement something sensible for this.
    self.assertTrue(is_chiral_implicit_hydrogen(c.top_back()));

    # All smiles variants must be identical
    usmi = m.unique_smiles()

    for i in range(10):
      m2 = LillyMolFromSmiles(m.random_smiles())
      self.assertEqual(m2.unique_smiles(), usmi)

  def test_invert_chirality(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("C[C@H](N)F"))
    usmi = m.unique_smiles()
    m.invert_chirality_on_atom(1);
    self.assertNotEqual(usmi, m.unique_smiles())
    m.remove_chiral_centre_at_atom(1);
    self.assertNotEqual(usmi, m.unique_smiles())
    self.assertFalse('@' in m.unique_smiles())

  def test_iterate_chiral_centres(self):
    m = Molecule();
    self.assertTrue(m.build_from_smiles("O[C@H]1[C@@H](O)C[C@@H](N)[C@H]1O CHEMBL268037"))
    atoms = Set_of_Atoms()
    for c in m.chiral_centres():
      atoms.append(c.atom())
    self.assertEqual(atoms, [1, 2, 5, 7])

  def test_charge(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC(=O)[O-]"))
    self.assertEqual(m.net_formal_charge(), -1)
    self.assertEqual(m.number_formal_charges(), 1)
    self.assertTrue(m.has_formal_charges())

    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC(=O)O"))
    self.assertEqual(m.formal_charge(3), 0)
    self.assertEqual(m.net_formal_charge(), 0)
    m.set_formal_charge(3, -1)
    self.assertEqual(m.formal_charge(3), -1)
    self.assertEqual(m.net_formal_charge(), -1)

  def test_xlogp(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O aspirin"))
    self.assertAlmostEqual(xlogp(m), 1.426)

  def test_alogp(self):
    m = Molecule()
    alogp = ALogP()
    self.assertTrue(m.build_from_smiles("C(C)[C@@H](NC1=NC=C2C(=C1)C(=C(C#N)C=N2)NC1=CC(=C(F)C=C1)Cl)C1=CC=CC=C1 CHEMBL197382"))
    self.assertAlmostEqual(alogp.logp(m), 6.601, places=3)

  def test_rotbond(self):
    m = Molecule()
    self.assertTrue(m.build_from_smiles("CC"))
    rotbond_calc = RotatableBonds()
    rotbond_calc.set_calculation_type(EXPENSIVE)
    self.assertEqual(rotbond_calc.rotatable_bonds(m), 0)
    self.assertTrue(m.build_from_smiles("CCC"))
    self.assertEqual(rotbond_calc.rotatable_bonds(m), 0)
    self.assertTrue(m.build_from_smiles("CC(F)(F)F"))
    self.assertEqual(rotbond_calc.rotatable_bonds(m), 0)
    self.assertTrue(m.build_from_smiles("CCCC"))
    self.assertEqual(rotbond_calc.rotatable_bonds(m), 1)
    self.assertTrue(m.build_from_smiles("C1CC1C"))
    self.assertEqual(rotbond_calc.rotatable_bonds(m), 0)

if __name__ == '__main__':
  absltest.main()

# Tests for set_of_atoms

import itertools
import random

from absl import app
from absl import logging
from absl.testing import absltest

from lillymol import *

class TestSetOfAtoms(absltest.TestCase):
  def test_empty(self):
    a = Set_of_Atoms()
    self.assertLen(a, 0)
    self.assertTrue(a.empty())
    self.assertEqual(a.size(), 0)

  def test_single_atom(self):
    a = Set_of_Atoms([2]);
    self.assertLen(a, 1)
    self.assertEqual(a[0], 2)
    self.assertEqual(a.size(), 1)

  def test_from_vector(self):
    for i in range(1, 10):
      v = list(range(0, i))
      a = Set_of_Atoms(v)
      self.assertEqual(a, v)

  def test_any_atoms_in_common(self):
    s1 = Set_of_Atoms([1, 2])
    self.assertTrue(s1.any_atoms_in_common(s1))
    self.assertEqual(s1.atoms_in_common(s1), 2)
    s2 = Set_of_Atoms([3])
    self.assertTrue(s2.any_atoms_in_common(s2))
    self.assertEqual(s2.atoms_in_common(s2), 1)
    self.assertFalse(s1.any_atoms_in_common(s2))
    self.assertFalse(s2.any_atoms_in_common(s1))

    s2 = Set_of_Atoms([1])
    self.assertTrue(s1.any_atoms_in_common(s2))
    self.assertTrue(s2.any_atoms_in_common(s1))
    self.assertEqual(s1.first_atom_in_common(s2), 1)
    self.assertEqual(s2.first_atom_in_common(s1), 1)
    self.assertEqual(s1.atoms_in_common(s2), 1)
    self.assertEqual(s2.atoms_in_common(s1), 1)

    s1 = Set_of_Atoms([0, 1, 2, 3])
    s2 = Set_of_Atoms([2, 3, 4, 5])
    self.assertTrue(s1.any_atoms_in_common(s2))
    self.assertTrue(s2.any_atoms_in_common(s1))
    self.assertEqual(s1.atoms_in_common(s2), 2)
    self.assertEqual(s2.atoms_in_common(s1), 2)
    self.assertEqual(s1.first_atom_in_common(s2), 2)
    self.assertEqual(s2.first_atom_in_common(s1), 2)

    # acknowledge that random in a unit test is never a good idea.
    random.shuffle(s2)
    self.assertEqual(s1.first_atom_in_common(s2), 2)

  def test_in(self):
    s = Set_of_Atoms(range(0, 10))
    for i in range(0, 10):
      self.assertTrue(i in s)

  def test_getitem(self):
    s = Set_of_Atoms(range(0, 10))
    for i in range(0, 10):
      self.assertEqual(s[i], i)

  def test_assignment(self):
    a = Set_of_Atoms([0, 1])
    a[0] = 5
    self.assertEqual(a, [5, 1])
    a[1] = 5
    self.assertEqual(a, [5, 5])

  def test_contains_both(self):
    a = Set_of_Atoms(range(0, 10))
    for (i, j) in itertools.combinations(range(0, 10), 2):
      self.assertTrue(a.contains_both(i, j));

  def test_extend(self):
    s1 = Set_of_Atoms(range(0, 5))
    s1.extend(range(5, 10))
    self.assertLen(s1, 10)
    for i in range(0, 10):
      self.assertEqual(s1[i], i)

  def test_iterate(self):
    s = Set_of_Atoms(range(0, 10))
    ndx = 0
    for atom in s:
      self.assertEqual(atom, ndx)
      ndx += 1

  def test_enumerate(self):
    s1 = Set_of_Atoms(range(0, 10))
    for ndx, atom in enumerate(s1):
      self.assertEqual(atom, ndx)

  def test_plus_equals(self):
    s1 = Set_of_Atoms(range(0, 5))
    s2 = Set_of_Atoms(range(5, 10))
    s1 += s2
    for ndx, atom in enumerate(s1):
      self.assertEqual(atom, ndx)


if __name__ == '__main__':
  absltest.main()

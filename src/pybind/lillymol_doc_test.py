"""Check that the python documentation describes methods that actually exist.

The API reference in docs/python is a set of markdown tables, one per class. A
table row naming a method that has since been renamed is worse than no
documentation at all - it sends the reader confidently down a dead end, and
nothing about it looks wrong on the page. This test reads those tables and asserts
every method named in them is really bound.

It scans every .md file in docs/python and keys off the table headings, so it
keeps working if the reference is split into other files, or if tables are added
for classes not listed here.

It deliberately does NOT assert the converse. Requiring every binding to be
documented would fail on every new method, which would train people to ignore it.
Undocumented methods are reported at the end as a to-do list instead.
"""

import os
import re

from absl import logging
from absl.testing import absltest

from lillymol import Atom, Bond, Molecule, Ring, Set_of_Atoms

# Table heading -> the class it describes. A heading not listed here is skipped,
# so an unrelated table cannot fail the test.
HEADING_TO_CLASS = {
    'Molecule Methods': Molecule,
    'Atom Methods': Atom,
    'Bond Methods': Bond,
    'Set_of_Atoms Methods': Set_of_Atoms,
    'Ring Methods': Ring,
}

# Rows that are not method names: the header, the rule, and the operator rows
# which document syntax rather than an attribute.
_NOT_A_METHOD = re.compile(r'^(Method|-+|m1\s|\\?_)')

# Leading identifier of a table cell, e.g. 'natoms(atomic_number)' -> 'natoms'.
_IDENTIFIER = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)')

_HEADING = re.compile(r'^#{2,3}\s+(.+?)\s*$')


def docs_directory():
  """The docs/python belonging to the same checkout as this test.

  Deliberately derived from __file__ rather than LILLYMOL_HOME. If you have more
  than one checkout, LILLYMOL_HOME may well point at a different one, and a test
  that validates someone else's documentation against this tree's bindings is
  worse than useless - it fails for reasons you cannot see in your own diff.
  """
  # src/pybind/lillymol_doc_test.py -> ../../docs/python
  here = os.path.dirname(os.path.abspath(__file__))
  return os.path.normpath(os.path.join(here, '..', '..', 'docs', 'python'))


def documented_methods(directory):
  """Yield (file, heading, cell, method_name) for every method table row."""
  for fname in sorted(os.listdir(directory)):
    if not fname.endswith('.md'):
      continue
    path = os.path.join(directory, fname)
    with open(path) as inp:
      heading = None
      for line in inp:
        matched = _HEADING.match(line)
        if matched:
          heading = matched.group(1)
          continue
        if heading not in HEADING_TO_CLASS or not line.startswith('|'):
          continue
        cell = line.split('|')[1].strip()
        if not cell or _NOT_A_METHOD.match(cell):
          continue
        identifier = _IDENTIFIER.match(cell)
        if identifier:
          yield fname, heading, cell, identifier.group(1)


class TestDocumentedMethodsExist(absltest.TestCase):

  def test_every_documented_method_is_bound(self):
    directory = docs_directory()
    self.assertTrue(os.path.isdir(directory), f'no docs directory {directory}')

    checked = 0
    wrong = []
    for fname, heading, cell, method in documented_methods(directory):
      checked += 1
      cls = HEADING_TO_CLASS[heading]
      if not hasattr(cls, method):
        # Offer the nearest bound names, which is almost always the rename.
        similar = sorted(a for a in dir(cls)
                         if not a.startswith('_') and a[:5] == method[:5])
        wrong.append(f'  {fname}, "{heading}": documented "{cell}" but '
                     f'{cls.__name__} has no "{method}". '
                     f'Similar: {similar if similar else "nothing"}')

    self.assertGreater(checked, 100,
                       'found almost no documented methods - has the reference '
                       'moved, or the table format changed?')
    self.assertEqual(wrong, [], 'documented methods that are not bound:\n' +
                     '\n'.join(wrong))
    logging.info('checked %d documented methods', checked)

  def test_report_undocumented_methods(self):
    """Not an assertion. Prints what is bound but never mentioned."""
    directory = docs_directory()
    mentioned = set()
    for fname in os.listdir(directory):
      if not fname.endswith('.md'):
        continue
      with open(os.path.join(directory, fname)) as inp:
        text = inp.read()
      mentioned |= set(re.findall(r'([A-Za-z_][A-Za-z0-9_]*)\s*\(', text))

    for cls in (Molecule, Atom, Bond, Set_of_Atoms, Ring):
      bound = {a for a in dir(cls) if not a.startswith('_')}
      undocumented = sorted(bound - mentioned)
      logging.info('%s: %d of %d methods are not mentioned in docs/python%s',
                   cls.__name__, len(undocumented), len(bound),
                   ': ' + ', '.join(undocumented) if undocumented else '')


if __name__ == '__main__':
  absltest.main()

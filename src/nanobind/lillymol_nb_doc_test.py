"""Check that docs/python documents nanobind-bound methods that exist.

The API reference in docs/python is a set of markdown tables, one per class. A
row naming a method that has since been renamed is worse than no documentation,
because it sends the reader confidently down a dead end.

This test deliberately does not require every bound method to be documented. It
only verifies that documented methods exist on the nanobind classes. Undocumented
methods are logged as a to-do list.
"""

import os
import re
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))

import lillymol

HEADING_TO_CLASS = {
    'Molecule Methods': lillymol.Molecule,
    'Atom Methods': lillymol.Atom,
    'Bond Methods': lillymol.Bond,
    'Set_of_Atoms Methods': lillymol.Set_of_Atoms,
    'Ring Methods': lillymol.Ring,
}

_NOT_A_METHOD = re.compile(r'^(Method|-+|m1\s|\\?_)')
_IDENTIFIER = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)')
_HEADING = re.compile(r'^#{2,3}\s+(.+?)\s*$')


def docs_directory():
  """Return docs/python for this checkout.

  Direct execution from src/nanobind can find ../../docs/python. Bazel tests run
  from runfiles and cannot include files outside the src workspace, so they need
  LILLYMOL_HOME to point at the repository root.
  """
  here = os.path.dirname(os.path.abspath(__file__))
  real_here = os.path.dirname(os.path.realpath(__file__))
  candidates = [
      os.path.normpath(os.path.join(here, '..', '..', 'docs', 'python')),
      os.path.normpath(os.path.join(real_here, '..', '..', 'docs', 'python')),
  ]
  home = os.environ.get('LILLYMOL_HOME')
  if home:
    candidates.append(os.path.join(home, 'docs', 'python'))
  else:
    raise FileNotFoundError("No LILLYMOL_HOME")

  print(candidates)
  for candidate in candidates:
    print(f"Checking documentation in {candidate}")
    if os.path.isdir(candidate):
      return candidate

  raise FileNotFoundError(f'Cannot find docs/python; set LILLYMOL_HOME ${LILLYMOL_HOME} to the repository root')


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


class TestNanobindDocumentedMethodsExist(unittest.TestCase):

  def test_every_documented_method_is_bound(self):
    directory = docs_directory()
    checked = 0
    wrong = []
    for fname, heading, cell, method in documented_methods(directory):
      checked += 1
      cls = HEADING_TO_CLASS[heading]
      if not hasattr(cls, method):
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
    print(f'checked {checked} documented methods in {directory}')

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

    for cls in (lillymol.Molecule, lillymol.Atom, lillymol.Bond,
                lillymol.Set_of_Atoms, lillymol.Ring):
      bound = {a for a in dir(cls) if not a.startswith('_')}
      undocumented = sorted(bound - mentioned)
      suffix = ': ' + ', '.join(undocumented) if undocumented else ''
      print(f'{cls.__name__}: {len(undocumented)} of {len(bound)} methods are not mentioned in docs/python{suffix}')


if __name__ == '__main__':
  unittest.main()

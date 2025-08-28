# Compare two smiles files, the first one contains multiple fragments.
# The second one contains the largest fragment.
# Identify the fragments and write a possible saltfile

import sys

from absl import app
from absl import logging

from lillymol import *
from lillymol_standardise import *

chemical_standardisation = Standardise()

def usage(rc :int):
  sys.exit(rc)

# smi1 contains fragments, smi2 is a single smiles of the parent.
def fragments(found, smi1, smi2):
  # Convert all to unique smiles
  m = MolFromSmiles(smi2)
  m.remove_isotopes()
  m.remove_all_chiral_centres()
  chemical_standardisation.process(m)
  parent = m.unique_smiles()

  for s in smi1.split('.'):
    # try a string match, might work.
    if s == smi2:
      continue

    m = MolFromSmiles(s)
    m.remove_isotopes()
    m.remove_all_chiral_centres()
    chemical_standardisation.process(m)
    if m.unique_smiles() == parent:  
      continue

    if m.natoms() > 20:
      continue

    usmi = m.unique_smiles()
    if usmi in found:
      found[usmi] += 1
    else:
      found[usmi] = 1

def counterions(found, inp1, inp2):
  for line1 in inp1:
    line2 = inp2.readline()
    fragments(found, line1.rstrip(), line2.rstrip())

def identify_counterions(argv):

  if len(argv) != 3:
    logging.error("Must specify exactly two input files")
    usage(1)

  found = {}

  set_auto_create_new_elements(True)
  chemical_standardisation.activate_all()

  with open(argv[1], 'r') as inp1:
    with open(argv[2], 'r') as inp2:
      counterions(found, inp1, inp2)

  for k, v in found.items():
    mol = MolFromSmiles(k)
    # Skip large and rare fragments
    if v < 5 and mol.natoms() > 10:
      continue
    print(f"{k} {mol.natoms()} {v}")

if __name__ == '__main__':
  app.run(identify_counterions)

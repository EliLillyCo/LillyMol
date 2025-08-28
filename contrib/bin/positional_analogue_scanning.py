#!/usr/bin/envy python

import sys
from itertools import combinations

from absl import app
from absl import flags
from absl import logging

FLAGS = flags.FLAGS
flags.DEFINE_string("smarts", "[cH]", "smarts for attachment point")
flags.DEFINE_integer("num_sub", 1, "Number of simultaneous additions")

from lillymol import *
from lillymol_io import *
from lillymol_query import *

def usage():
  sys.exit(0)

# Implementation ideas copied from Pat Walters
# https://practicalcheminformatics.blogspot.com/2020/04/positional-analogue-scanning.html#:~:text=By%20generating%20a%20set%20of%20positional%20analogs%2C%20one,where%20PAS%20can%20be%20combined%20with%20computational%20analysis.

duplicates_skipped = 0
def attach_atom(mol_in:Molecule, atomic_number:int, query:SubstructureQuery, num_sub:int, seen:set):
  """Attach an atom with atomic number `atomic_number` to each atom in `mol` that
  matches `query`
    seen:       a set of unique smiles used for uniqueness determinations
    Returns: A list of the newly formed molecules.
  """

  result = []    # to be returned

  sresults = SubstructureResults()
  if query.substructure_search(mol_in, sresults) == 0:
    print("No query match\n", file=sys.stderr)
    return result

  match_atoms = [x[0] for x in sresults]
  n_combos = combinations(match_atoms, num_sub)
  print(n_combos,file=sys.stderr)
  seq = 0
  for combo in n_combos:
    newmol = Molecule(mol_in)
    print(f"{mol_in.name()} {combo}",file=sys.stderr)
    for idx in combo:
      newmol.add_atom(atomic_number)
      newmol.add_bond(idx, newmol.natoms() - 1, SINGLE_BOND)

    usmi = newmol.unique_smiles()
    if usmi in seen:
      duplicates_skipped += 1
      continue

    seen.add(usmi)
    newmol.set_name(f"{mol_in.name()}.{seq}")
    result.append(newmol)
    seq += 1

  return result


def main(argv):
  if len(argv) < 2:
    logging.warning("Must specify an input file")
    usage()

  query = SubstructureQuery()
  if not query.build_from_smarts(FLAGS.smarts):
    logging.error("Invalid smarts %s", FLAGS.smarts)

  seen = set()
  with ReaderContext(argv[1]) as reader:
    for mol in reader:
      print(mol)
      newmols = attach_atom(mol, 9, query, FLAGS.num_sub, seen)
      for x in newmols:
        print(x)

  print(f"Skipped {duplicates_skipped} duplicates")


if __name__ == '__main__':
  app.run(main)


#!/usr/bin/env python

# Demonstrates descriptor computation via the IWDescr object.

from absl import app
from absl import logging

from lillymol_io import ReaderContext
from lillymol_tools import IWDescr


# Write smiles, id, AMW and longest path for `mol`.
# name_to_col is a mapping from feature names to column numbers.
# See https://github.com/IanAWatson/LillyMol/blob/master/docs/Molecule_Tools/iwdescr.md
# for feature definitions.
def compute_descriptors(iwdescr, name_to_col, mol):
  descriptors = iwdescr.process(mol)
  amw = name_to_col["amw"]
  mxdst = name_to_col["mxdst"]
  print(f"{mol.smiles()} {mol.name()} {descriptors[amw]:.2f} {int(descriptors[mxdst])}")

def main(argv):
  if len(argv) == 1:
    logging.error("Must specify input file")
    return

  iwdescr = IWDescr()

  name_to_col = {
      name: col for col, name in enumerate(iwdescr.feature_names())
  }

  with ReaderContext(argv[1]) as input:
    for mol in input:
      compute_descriptors(iwdescr, name_to_col, mol)

if __name__ == "__main__":
  app.run(main)

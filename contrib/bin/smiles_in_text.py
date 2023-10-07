# Scan through files and extract all tokens interpretable as smiles.
# This is part of the open source release of LillyMol to ensure there
# is no proprietary structural information disclosed.
# Files are scanned for strings that can be interpreted as smiles
# and this tool reports valid structures/smiles.
# Those structures are then looked up in Pubchem, Chembl and other
# public data sources. Then those that remain are looked up in
# the Lilly database and we find no structures.

# in_pubchem.sh -i smi -i ICTE -U - structures.smi | \
# in_chembl.sh -U - -i smi -v - | \
# ... | \
# in_lilly_database.sh -F - -i smi -

# nothing is found.

import os
import re
import string

from absl import app

from lillymol import *

def main(argv):
  """Scan through all files below the current directory and extract smiles"""
  # https://stackoverflow.com/questions/18429143/strip-punctuation-with-regex-python
  # p = re.compile("[" + re.escape(string.punctuation) + "]")
  my_punct = ['!', '"', '$', '&', "'", '*', ',', 
              ';', '<', '>', '?', '^', '_', 
              '`', '{', '|', '}', '~', '“', '”']
  p = re.compile("[" + re.escape("".join(my_punct)) + "]")
  do_not_process = re.compile(".*\.(png|pdf)$")

  # initially implemented this to avoid finding the same molecule
  # multiple times, but we actually want to find all instances of
  # problematic structures.
  seen = set()

  for root, dirs, files in os.walk("."):
    for name in files:
      fname = os.path.join(root, name)
      if do_not_process.match(fname):
        continue
      with open(fname, "r") as reader:
        for line in reader:
          no_punctuation = p.sub(" ", line.rstrip())
          for token in no_punctuation.split():
            # Ignore small 'molecules'
            if len(token) < 10:
              continue
            mol = MolFromSmiles(token)
            if not mol:
              continue
            if mol.natoms() < 10:
              continue
            # isotopes are used for tests and documentation
            if mol.number_isotopic_atoms() > 0:
              continue
            # do not update `seen` because we want all instances
            # if mol.unique_smiles() in seen:
            #   continue
            # seen.add(mol.unique_smiles())
            print(f"{mol.aromatic_smiles()} {fname}")

if __name__ == "__main__":
  app.run(main)

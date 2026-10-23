# Compute RDKit 2D descriptors
# Takes significantly longer than iwdescr.

import numpy as np

from absl import app
from absl import logging

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

def compute_descriptors(mol, name, calc):
  """Compute RDKit 2D descriptors for `mol` and write.
   Args:
     mol: molecule
     name: name of the molecule
     calc: descriptors to be computed
  """
  ds = calc.CalcDescriptors(mol)
  print(name, "", end="")
  # Must be a better way of doing this output.
  # would also like to not print so many significant digits.
  x = [str(v) for v in ds]
  print(' '.join(x))


def rdkit_descriptors(argv):
  if len(argv) == 1:
    logging.error("Must specify smiles file as an argument")
    return 1

  calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
  names = calc.GetDescriptorNames()
  print("ID", " ".join(names))
  with open(argv[1], "r") as input:
    for line in input:
      f = line.rstrip().split(' ')
      mol = Chem.MolFromSmiles(f[0])
      if mol is None:
        logging.error("Bad smiles %s", f[0])
        continue

      compute_descriptors(mol, f[1], calc)


if __name__ == "__main__":
  app.run(rdkit_descriptors)

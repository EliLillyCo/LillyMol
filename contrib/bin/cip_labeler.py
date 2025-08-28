# Label stereo centres with rdCIPLabeler
# R and S sites are assigned an isotopic label according to R or S.

from dataclasses import dataclass
import sys

from absl import app
from absl import flags
from absl import logging


FLAGS = flags.FLAGS

flags.DEFINE_integer("R", 8, "Isotope for R centres")
flags.DEFINE_integer("S", 9, "Isotope for R centres")
flags.DEFINE_boolean("f", False, "Work as a TDT filter - used internally by gfp_make")
flags.DEFINE_boolean("descriptors", False, "Generate a descriptor file")
flags.DEFINE_boolean("verbose", False, "Verbose output")

from rdkit import Chem

@dataclass
class Options:
  cip_r: int = 8
  cip_s: int = 9
  molecules_read: int = 0
  molecules_with_chiral_centres: int = 0

def usage(rc):
  print("Applies isotopic labels according to R or S stereo configurations")
  print("Generates descriptors of R and S counts with the --descriptors option")
  print("Use --help for more info")
   
  sys.exit(rc)

def cip_label(options:Options, smiles:str):
  # Return a smiles with isotopic labels on stereo sites in `smiles`.
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    logging.error("Invalid smiles %s", smiles)
    return smiles

  options.molecules_read += 1

  Chem.rdCIPLabeler.AssignCIPLabels(mol)
  stereo_this_molecule = 0
  for atom in mol.GetAtoms():
    c = atom.GetChiralTag()
    if not atom.HasProp("_CIPCode"):
      continue

    stereo_this_molecule += 1
    cipcode = atom.GetProp("_CIPCode")
    if cipcode == 'S':
      atom.SetIsotope(options.cip_s)
    elif cipcode == 'R':
      atom.SetIsotope(options.cip_r)
    
  if stereo_this_molecule == 0:
    return smiles

  options.molecules_with_chiral_centres += 1
  return Chem.MolToSmiles(mol)

def cip_label_filter(options:Options, line:str):
  """Reading a fingerprint file. Find the $SMI<> records and process them
  """
  if not line.startswith("$SMI<"):
    print(line)
    return

  smi = line[5:-1]

  smi = cip_label(options, smi)
  print(f"$SMI<{smi}>")

def cip_label_smiles_line(options:Options, smiles_id:str):
  """Write isotopically labelled smiles for the smiles in `smiles_id`.
    Args:
      options:
      smiles_id: sring containing smiles and id.
  """
  f = smiles_id.split()
  smiles = cip_label(options, f[0])
  print(f"{smiles} {' '.join(f[1:])}")

def cip_label_descriptors(options:Options, line:str):
  """Generate descriptors for the smiles in `line` and write.
  """
  f = line.split()
  mol = Chem.MolFromSmiles(f[0])
  if mol is None:
    logging.error("Invalid smiles %s", line)
    print(f"{f[1]} 0 0")
    return

  options.molecules_read += 1

  Chem.rdCIPLabeler.AssignCIPLabels(mol)
  sumr = 0
  sums = 0
  for atom in mol.GetAtoms():
    c = atom.GetChiralTag()
    if not atom.HasProp("_CIPCode"):
      continue

    cipcode = atom.GetProp("_CIPCode")
    if cipcode == 'S':
      sums += 1
    elif cipcode == 'R':
      sumr += 1

  if sumr > 0 or sums > 0:
    options.molecules_with_chiral_centres += 1

  print(f"{f[1]} {sumr} {sums}")


def main(argv):
  if len(argv) == 1:
    logging.error("Must provide input file")
    usage(1)

  options = Options()

  options.cip_r = FLAGS.R
  options.cip_s = FLAGS.S

  if FLAGS.descriptors:
    print("Name cip_R cip_S")
    for fname in argv[1:]:
      with open(fname, 'r') as input:
        for line in input:
          cip_label_descriptors(options, line.rstrip())
  elif FLAGS.f:
    for line in sys.stdin:
      cip_label_filter(options, line.rstrip())

  elif argv[1] == '-':
    for line in sys.stdin:
      cip_label_smiles_line(options, line.rstrip())

  else:
    for fname in argv[1:]:
      with open(fname, 'r') as input:
        for line in input:
          cip_label_smiles_line(options, line.rstrip())

  if FLAGS.verbose:
    logging.info("Read %d molecules", options.molecules_read)
    logging.info("%d had stereo centres", options.molecules_with_chiral_centres);

if __name__ == "__main__":
  app.run(main)

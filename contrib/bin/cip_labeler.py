#!/usr/bin/env python
"""Mark centres using their original R/S CIP assignments.

Default: apply isotope 8 to R centres and isotope 9 to S centres.
Element mappings take precedence over isotope labelling for matching centres:
  cip_labeler.py --R_element=C:Cr --S_element=C:Cs input.smi
Repeat either mapping option to handle additional source elements. Unmapped
centres retain isotope labelling. Replacements preserve existing isotopes,
bonds, formal charges and hydrogen counts; descriptors count original labels.
The replacement elements are structural markers, not chemically validated
products. For example, tetra-coordinate Cs will fail RDKit valence sanitisation;
RDKit consumers of such output must read with sanitize=False.
"""

from dataclasses import dataclass, field
import sys

from absl import app
from absl import flags
from absl import logging


FLAGS = flags.FLAGS

flags.DEFINE_integer("R", 8, "Isotope for R centres")
flags.DEFINE_integer("S", 9, "Isotope for S centres")
flags.DEFINE_multi_string("R_element", [], "Replace R centres: SOURCE:TARGET, e.g. C:Cr")
flags.DEFINE_multi_string("S_element", [], "Replace S centres: SOURCE:TARGET, e.g. C:Cs")
flags.DEFINE_boolean("f", False, "Work as a TDT filter - used internally by gfp_make")
flags.DEFINE_boolean("descriptors", False, "Generate a descriptor file")
flags.DEFINE_boolean("verbose", False, "Verbose output")

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

@dataclass
class Options:
  cip_r: int = 8
  cip_s: int = 9
  r_elements: dict[int, int] = field(default_factory=dict)
  s_elements: dict[int, int] = field(default_factory=dict)
  molecules_read: int = 0
  molecules_with_chiral_centres: int = 0

def usage(rc):
  print("Applies isotopic labels or element replacements according to R/S configurations")
  print("Generates descriptors of R and S counts with the --descriptors option")
  print("Use --help for more info")
   
  sys.exit(rc)

def parse_element_mappings(specifications):
  """Validate SOURCE:TARGET symbols and return atomic-number mappings."""
  periodic_table = Chem.GetPeriodicTable()
  elements = {periodic_table.GetElementSymbol(n): n for n in range(1, 119)}
  mappings = {}
  for specification in specifications:
    fields = specification.split(':')
    if len(fields) != 2 or not all(fields):
      raise ValueError(f"Invalid element mapping {specification!r}; expected SOURCE:TARGET")
    numbers = []
    for symbol in fields:
      number = elements.get(symbol)
      if number is None:
        raise ValueError(f"Invalid element symbol {symbol!r}")
      numbers.append(number)
    source, target = numbers
    if source in mappings and mappings[source] != target:
      raise ValueError(f"Conflicting replacements for {fields[0]}")
    mappings[source] = target
  return mappings


def find_cip_centres(mol):
  """Snapshot original atom indices, elements and uppercase R/S assignments."""
  rdCIPLabeler.AssignCIPLabels(mol)
  centres = []
  for atom in mol.GetAtoms():
    if not atom.HasProp("_CIPCode"):
      continue
    label = atom.GetProp("_CIPCode")
    if label in ('R', 'S'):
      centres.append((atom.GetIdx(), atom.GetAtomicNum(), label))
  return centres


def cip_label(options:Options, smiles:str):
  """Identify original CIP centres before applying any changes to the molecule."""
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    logging.error("Invalid smiles %s", smiles)
    return smiles

  options.molecules_read += 1
  centres = find_cip_centres(mol)
  if not centres:
    return smiles

  options.molecules_with_chiral_centres += 1
  # No CIP reassignment during this pass: replacements can change priorities at
  # other centres. All decisions use the snapshot of the input molecule.
  for index, original_element, label in centres:
    atom = mol.GetAtomWithIdx(index)
    replacements = options.r_elements if label == 'R' else options.s_elements
    replacement = replacements.get(original_element)
    if replacement is not None:
      # Preserve attached hydrogens explicitly rather than asking the new
      # element's valence rules to infer a different hydrogen count.
      hydrogens = atom.GetTotalNumHs()
      atom.SetAtomicNum(replacement)
      atom.SetNumExplicitHs(hydrogens)
      atom.SetNoImplicit(True)
    else:
      atom.SetIsotope(options.cip_r if label == 'R' else options.cip_s)

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

  centres = find_cip_centres(mol)
  sumr = sum(label == 'R' for _, _, label in centres)
  sums = sum(label == 'S' for _, _, label in centres)

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
  try:
    options.r_elements = parse_element_mappings(FLAGS.R_element)
    options.s_elements = parse_element_mappings(FLAGS.S_element)
  except ValueError as error:
    raise app.UsageError(str(error)) from error


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

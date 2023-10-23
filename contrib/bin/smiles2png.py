""" Given a smiles file, create a set of .png files which are sent to eog."""

from dataclasses import dataclass

import re
import os
import sys
from typing import List, Optional

from absl import app
from absl import flags
from absl import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

FLAGS = flags.FLAGS

flags.DEFINE_multi_string("input", None, "Smiles file, or use argv")
flags.DEFINE_multi_string("smiles", None, "Individual smiles to depict")
flags.DEFINE_multi_string("align", None, "Smiles of structure to align")
flags.DEFINE_string("stem", None, "File name stem for .png files created")
flags.DEFINE_integer("n", 30, "Number of smiles to process")
flags.DEFINE_integer("mols_per_row", 1, "Generate a grid plot with this many molecules per row")
flags.DEFINE_integer("nrows", 1, "Number of rows on each page")
flags.DEFINE_integer("x", 300, "Image size X")
flags.DEFINE_integer("y", 300, "Image size Y")
flags.DEFINE_boolean("keep", False, "Keep the underlying .png file(s)")
flags.DEFINE_boolean("verbose", False, "Verbose output")

UNAME = "_Name"

@dataclass
class Smiles2PngConfig:
  # pylint: disable=too-many-instance-attributes
  """Configurable parameters for plotting."""
  nplot: int = 30
  mols_per_row: int = 1
  rows_per_page: int = 1
  plot_x: int = 300
  plot_y: int = 300
  align: Optional[List[Chem.RWMol]] = None
  keep_png: bool = False
  verbose: bool = False

# pylint: disable=line-too-long

def generate_aligned_coordinates(mols: List[Chem.rdchem.Mol],
                                 align: List):
  """Generate aligned 2D coordinates for the molecules in `mols`.
  Args:
  Returns:
  """
  no_match = 0
  for mol in mols:
    # might be a clever way to do this with an any() construct...
    got_match = False
    for a in align: # pylint: disable=invalid-name
      if mol.HasSubstructMatch(a):
        AllChem.GenerateDepictionMatching2DStructure(mol, a)
        got_match = True
        break
    if not got_match:
      AllChem.Compute2DCoords(mol)
      no_match += 1

  if no_match > 0:
    print(f"Warning {no_match} of {len(mols)} did not match the alignment template", file=sys.stderr)

def generate_plots_grid(mols: List[Chem.rdchem.Mol],
                        png_stem: str,
                        config: Smiles2PngConfig) -> List[str]:
  """Convert molecules in `mols` to png files.
  this version creates a grid plot, according to the settings in `config`.
  Args:
    mols: a list of RDKit molecules
    png_stem: a file name stem for creating the .png files
    config: settings for the processing.
  Returns:
    A list of the .png files created
  """
  ndx = 0
  pngs = []
  per_page = config.mols_per_row * config.rows_per_page
  for i in range(0, len(mols), per_page):
    gridmols = mols[i : i + per_page]
    ids = [m.GetProp(UNAME) for m in gridmols]
    print(f"Drawing {len(gridmols)} molecules, per for {config.mols_per_row}")
    png = Draw.MolsToGridImage(gridmols, molsPerRow=config.mols_per_row, legends=ids, returnPNG=True)
    png_fname = f"{png_stem}{ndx}.png"
    with open(png_fname, "wb") as output:
      output.write(png)
    pngs.append(png_fname)
    ndx += 1

  return pngs

def generate_plots(mols: List[Chem.rdchem.Mol],
                   png_stem: str,
                   config: Smiles2PngConfig) -> List[str]:
  """Generate single molecule at a time plots of the molecules in `mols`.

  Args:
    mols:
    png_stem: used for generating .png files
    config:
  Returns:
    A list of the .png files produced.
  """
  pngs: List[str] = []  # to be returned.
  for (ndx, mol) in enumerate(mols):
    if len(mols) == 1:
      png_fname = f"{png_stem}.png"
    else:
      png_fname = f"{png_stem}{ndx}.png"
    if config.verbose:
      print(f"Drawing plot for {Chem.MolToSmiles(mol)}")
    Draw.MolToFile(mol, png_fname, size=(config.plot_x, config.plot_y), legend=mol.GetProp(UNAME))

    pngs.append(png_fname)

  return pngs

def display_pngs(pngs: List[str],
                 config: Smiles2PngConfig):
  """Given a list of png files, pass to eog for viewing.

  If specified in `config`, remove them when done
  Args:
    pngs: list of .png files to be passed to eog.
    config:
  """
  if len(pngs) == 0:
    print("No png files produced")
    return

  files = ' '.join(pngs)
  os.system(f"eog {files}")

  if not config.keep_png:
    for png in pngs:
      os.remove(png)


def read_smiles(input_stream,
                lines: List[str],
                config: Smiles2PngConfig) -> int:
  """Read lines from `input` and append to `lines`.

  Stop reading if/when config.nplot is reached.
  Args:
    input_stream:
    lines:
  Returns:
    length(lines)
  """
  for line in input_stream:
    lines.append(line.rstrip())
    if len(lines) >= config.nplot:
      return len(lines)

  return len(lines)

def get_smiles(smiles_fname: str,
               lines: List[str],
               config: Smiles2PngConfig) -> Optional[int]:
  """Wrapper for read_smiles depending on `smiles_fname`."""
  if smiles_fname == '-':
    return read_smiles(sys.stdin, lines, config)
  with open(smiles_fname, "r") as input_stream:
    return read_smiles(input_stream, lines, config)

def do_smiles2png(smiles_fname: List[str],
                  smiles_on_command_line: Optional[List[str]],
                  png_stem: str,
                  config: Smiles2PngConfig):
  """Generate .png files from the molecules in `smiles_fname`.
  Args:
    smiles_fname: list of files containing smiles
    png_stem: used for generating .png files
    config:
  """
  lines: List[str] = []
  if smiles_on_command_line is not None:
    lines = smiles_on_command_line
  else:
    for fname in smiles_fname:
      get_smiles(fname, lines, config)

  if config.verbose:
    print(f"Read {len(lines)} lines", file=sys.stderr)

  if len(lines) == 0:
    print("No smiles", file=sys.stderr)
    return

  mols = []
  for line in lines:
    f = line.split(' ') # pylint: disable=invalid-name
    mols.append(Chem.MolFromSmiles(f[0]))
    if len(f) > 1:
      mols[-1].SetProp(UNAME, f[1])
    else:
      mols[-1].SetProp(UNAME, "")

  if config.verbose:
    print(f"Generating plots for {len(mols)} molecules", file=sys.stderr)

  if config.align is not None:
    generate_aligned_coordinates(mols, config.align)
  else:
    _ = [AllChem.Compute2DCoords(m) for m in mols]

  if config.mols_per_row > 1:
    pngs = generate_plots_grid(mols, png_stem, config)
  else:
    pngs = generate_plots(mols, png_stem, config)


  if len(pngs) == 0:
    print("No png files produced")
    return

  display_pngs(pngs, config)

def generate_png_stem(fnames: List[str],
                      smiles_on_command_line: Optional[List[str]],
                      config: Smiles2PngConfig) -> str:
  """Generate a stem for generated png files.

  Complicated by the fact that input might be from stdin,
  and so the first file name might be '-'.

  Args:
    fnames: List of smiles files
    config:
  Returns:
    file name stem for png files.
  """
  if smiles_on_command_line or fnames[0] == '-':
    result =  f"/tmp/smi2png{os.getpid()}"
  else:
    result = re.sub(r'.smi$', "", fnames[0])

  if config.verbose:
    print(f"png stem {result}")

  return result


def smiles2png(argv):
  """Generate png from smiles using RDKit"""

  smiles_fname = []

  # Currently these are mutually exclusive...
  if FLAGS.input is not None:
    smiles_fname = FLAGS.input
  elif len(argv) > 0:
    smiles_fname = argv[1:]

  png_stem = FLAGS.stem

  config = Smiles2PngConfig()
  if FLAGS.align is not None:
    config.align = [Chem.MolFromSmiles(smi) for smi in FLAGS.align]
    for mol in config.align:
      AllChem.Compute2DCoords(mol)

  smiles_on_command_line: Optional[List[str]] = FLAGS.smiles

  config.nplot = FLAGS.n
  config.mols_per_row = FLAGS.mols_per_row
  config.rows_per_page = FLAGS.nrows
  config.plot_x = FLAGS.x
  config.plot_y = FLAGS.y
  config.keep_png = FLAGS.keep
  config.verbose = FLAGS.verbose

  if png_stem:
    config.keep_png = True
  else:
    png_stem = generate_png_stem(smiles_fname, smiles_on_command_line, config)

  do_smiles2png(smiles_fname, smiles_on_command_line, png_stem, config)

if __name__ == '__main__':
  app.run(smiles2png)

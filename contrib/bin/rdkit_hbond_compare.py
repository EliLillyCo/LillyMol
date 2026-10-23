# Compare LillyMol's RDKit compatible hydrogen bond counts against RDKit.
#
# LillyMol's Molecule::RDKitNumHAcceptors and RDKitNumHDonors are hand coded
# transcriptions of the SMARTS that RDKit compiles into CalcNumHBA and
# CalcNumHBD. There is no mechanism keeping the two in sync, so run this after
# an RDKit upgrade, or after touching the transcription, and look at the
# residual.
#
# Note that the comparison is against rdMolDescriptors.CalcNumHBA, NOT against
# the HAcceptorSmarts pattern in rdkit/Chem/Lipinski.py. Those two disagree -
# Lipinski.py says nH0 where the compiled descriptor says nH0X2 - and it is
# CalcNumHBA that NumHAcceptors actually calls.
#
# Some disagreement is expected and is not a bug. It is dominated by
# aromaticity perception, which LillyMol and RDKit will never fully share.
# The -a option reports how much of the residual is attributable to that.
#
#   rdkit_hbond_compare.py file.smi
#   rdkit_hbond_compare.py -a -f mismatches.smi file.smi

import sys

from absl import app
from absl import flags
from absl import logging

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolDescriptors

from lillymol import *

FLAGS = flags.FLAGS

flags.DEFINE_boolean("aromaticity", False,
                     "Also report how many mismatches involve a molecule where "
                     "the two toolkits disagree about aromatic atom count", short_name="a")
flags.DEFINE_string("failures", None,
                    "Write mismatching molecules to this smiles file", short_name="f")
flags.DEFINE_integer("report", 0,
                     "Report progress every this many molecules", short_name="r")


def aromatic_atom_count_rdkit(mol):
  return sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())


def aromatic_atom_count_lillymol(mol):
  return mol.aromatic_atom_count()


def compare(fname, failures):
  """Compare LillyMol and RDKit hydrogen bond counts over `fname`.

  Args:
    fname: smiles file
    failures: open file for mismatching molecules, or None
  Returns:
    True on success.
  """
  nmols = 0
  bad_smiles = 0
  agree_acceptor = 0
  agree_donor = 0
  agree_both = 0
  arom_differs_in_mismatch = 0
  arom_differs = 0

  with open(fname, "r") as inp:
    for line in inp:
      line = line.rstrip()
      if not line:
        continue
      tokens = line.split()
      smi = tokens[0]
      name = tokens[1] if len(tokens) > 1 else f"mol{nmols}"

      lmol = MolFromSmiles(smi)
      rmol = Chem.MolFromSmiles(smi)
      if lmol is None or rmol is None:
        bad_smiles += 1
        continue

      nmols += 1
      if FLAGS.report and nmols % FLAGS.report == 0:
        logging.info("Processed %d molecules", nmols)

      lhba = lmol.rdkit_num_h_acceptors()
      lhbd = lmol.rdkit_num_h_donors()
      rhba = rdMolDescriptors.CalcNumHBA(rmol)
      rhbd = rdMolDescriptors.CalcNumHBD(rmol)

      arom_same = True
      if FLAGS.aromaticity:
        arom_same = aromatic_atom_count_lillymol(lmol) == aromatic_atom_count_rdkit(rmol)
        if not arom_same:
          arom_differs += 1

      a_ok = lhba == rhba
      d_ok = lhbd == rhbd
      agree_acceptor += a_ok
      agree_donor += d_ok
      agree_both += a_ok and d_ok

      if a_ok and d_ok:
        continue

      if not arom_same:
        arom_differs_in_mismatch += 1
      if failures is not None:
        print(f"{smi} {name} lillymol {lhba} {lhbd} rdkit {rhba} {rhbd}", file=failures)

  if nmols == 0:
    logging.error("No molecules read from %s", fname)
    return False

  def pct(n):
    return f"{n:6d}  {100.0 * n / nmols:6.2f}%"

  print(f"{nmols} molecules, {bad_smiles} could not be parsed by both toolkits")
  print(f"  acceptors agree {pct(agree_acceptor)}")
  print(f"  donors agree    {pct(agree_donor)}")
  print(f"  both agree      {pct(agree_both)}")
  if FLAGS.aromaticity:
    mismatches = nmols - agree_both
    print(f"  aromatic atom count differs {pct(arom_differs)}")
    if mismatches:
      frac = 100.0 * arom_differs_in_mismatch / mismatches
      print(f"  of the {mismatches} mismatches, {arom_differs_in_mismatch} "
            f"({frac:.1f}%) are molecules where aromaticity also differs")

  return True


def main(argv):
  if len(argv) == 1:
    logging.error("Must specify a smiles file as an argument")
    return 1

  RDLogger.DisableLog("rdApp.*")

  failures = open(FLAGS.failures, "w") if FLAGS.failures else None
  try:
    for fname in argv[1:]:
      if not compare(fname, failures):
        return 1
  finally:
    if failures is not None:
      failures.close()

  return 0


if __name__ == "__main__":
  app.run(main)

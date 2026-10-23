#!/usr/bin/env python3
"""Compare LillyMol and RDKit atom hybridization assignments."""

import argparse
import sys

from lillymol import MolFromSmiles, hybridization

try:
  from rdkit import Chem
except ImportError as e:
  raise SystemExit("RDKit is required for this comparison script") from e


def tokens(line):
  fields = line.strip().split()
  if not fields:
    return None, None
  smiles = fields[0]
  name = fields[1] if len(fields) > 1 else smiles
  return smiles, name


def compare_record(smiles, name, verbose):
  lm = MolFromSmiles(smiles)
  rdkit = Chem.MolFromSmiles(smiles)
  if lm is None or rdkit is None:
    print(f"{name}: cannot parse {smiles}", file=sys.stderr)
    return 1

  if lm.natoms() != rdkit.GetNumAtoms():
    print(f"{name}: atom count mismatch LillyMol {lm.natoms()} RDKit {rdkit.GetNumAtoms()}")
    return 1

  diffs = 0
  atom_diffs = []
  for i in range(lm.natoms()):
    lhs = str(hybridization(lm, i))
    rhs = str(rdkit.GetAtomWithIdx(i).GetHybridization())
    if lhs.startswith("Hybridization."):
      lhs = lhs.split(".", 1)[1]
    if lhs != rhs:
      print(f"{name} atom {i} {lm.smarts_equivalent_for_atom(i)} LillyMol={lhs} RDKit={rhs}")
      diffs += 1
      atom_diffs.append(i)
    elif verbose:
      atom = rdkit.GetAtomWithIdx(i)
      print(f"{name} atom {i} {lm.smarts_equivalent_for_atom(i)} {lhs}")

  if diffs > 0:
    print(f"{lm.name()} {diffs} differences")
    for i in atom_diffs:
      lm.set_isotope(i, 1)
    print(f"{lm.smiles()} {name} {diffs} DIFFS")
  return diffs


def main(argv):
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("smiles", nargs="*", help="SMILES files. Reads stdin if omitted or '-' is present")
  parser.add_argument("-v", "--verbose", action="store_true", help="also report matching atoms")
  args = parser.parse_args(argv)

  inputs = args.smiles or ["-"]
  diffs = 0
  molecules_read = 0
  molecules_with_differences = 0
  for fname in inputs:
    stream = sys.stdin if fname == "-" else open(fname)
    with stream:
      for line in stream:
        molecules_read += 1
        smiles, name = tokens(line)
        if smiles is None or smiles.startswith("#"):
          continue
        tmp = compare_record(smiles, name, args.verbose)
        if tmp > 0:
          molecules_with_differences += 1
          diffs += tmp

  if diffs > 0:
    print(f"Read {molecules_read}, {molecules_with_differences} had differences, total {diffs} diffs found")

  return 1 if diffs else 0


if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))

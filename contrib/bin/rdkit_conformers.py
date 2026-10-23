#!/usr/bin/env python3
"""Generate 3D conformers from a whitespace-delimited SMILES file."""

import argparse
import sys
from pathlib import Path

try:
  from rdkit import Chem
  from rdkit.Chem import AllChem
except ImportError as exc:
  raise SystemExit("RDKit is required: install it before running this program") from exc


def parse_args():
  parser = argparse.ArgumentParser(
      description="Generate energy-ranked 3D conformers with RDKit ETKDGv3")
  parser.add_argument("input", help="SMILES file, or - for standard input")
  parser.add_argument("-o", "--output", required=True, help="Output SDF file")
  parser.add_argument("-n", "--num-conformers", type=int, default=20,
                      help="maximum conformers generated per molecule (default: 20)")
  parser.add_argument("--prune-rms", type=float, default=0.5,
                      help="RMS threshold in Angstrom for pruning similar conformers "
                           "(default: 0.5; negative disables pruning)")
  parser.add_argument("--seed", type=int, default=0xF00D,
                      help="random seed for reproducible embedding (default: 61453)")
  parser.add_argument("--threads", type=int, default=0,
                      help="threads used by RDKit; 0 uses all available threads")
  parser.add_argument("--max-iterations", type=int, default=1000,
                      help="maximum force-field iterations per conformer (default: 1000)")
  parser.add_argument(
      "--force-field", choices=("auto", "mmff94", "mmff94s", "uff", "none"),
      default="auto", help="optimization force field (default: auto: MMFF94s then UFF)")
  parser.add_argument("--remove-hydrogens", action="store_true",
                      help="remove explicit hydrogens from SDF output after optimization")
  args = parser.parse_args()

  if args.num_conformers < 1:
    parser.error("--num-conformers must be positive")
  if args.threads < 0:
    parser.error("--threads must be non-negative")
  if args.max_iterations < 1:
    parser.error("--max-iterations must be positive")
  return args


def read_smiles(fname):
  input_stream = sys.stdin if fname == "-" else open(fname, encoding="utf-8")
  try:
    for line_number, line in enumerate(input_stream, 1):
      line = line.strip()
      if not line or line.startswith("#"):
        continue
      fields = line.split(maxsplit=1)
      smiles = fields[0]
      name = fields[1] if len(fields) == 2 else f"molecule_{line_number}"
      yield line_number, smiles, name
  finally:
    if input_stream is not sys.stdin:
      input_stream.close()


def optimize(mol, args):
  """Returns (force field name, [(status, energy), ...])."""
  if args.force_field == "none":
    return "none", [(3, None)] * mol.GetNumConformers()

  requested = args.force_field
  if requested in ("auto", "mmff94", "mmff94s"):
    variant = "MMFF94" if requested == "mmff94" else "MMFF94s"
    if AllChem.MMFFHasAllMoleculeParams(mol):
      results = AllChem.MMFFOptimizeMoleculeConfs(
          mol, numThreads=args.threads, maxIters=args.max_iterations,
          mmffVariant=variant)
      return variant, results
    if requested != "auto":
      return variant, [(2, None)] * mol.GetNumConformers()

  if requested in ("auto", "uff"):
    if AllChem.UFFHasAllMoleculeParams(mol):
      results = AllChem.UFFOptimizeMoleculeConfs(
          mol, numThreads=args.threads, maxIters=args.max_iterations)
      return "UFF", results
    return "UFF", [(2, None)] * mol.GetNumConformers()

  raise AssertionError(f"Unhandled force field {requested}")


def conformers(mol, args):
  mol = Chem.AddHs(mol)
  params = AllChem.ETKDGv3()
  params.randomSeed = args.seed
  params.pruneRmsThresh = args.prune_rms
  params.numThreads = args.threads
  params.useRandomCoords = False

  conformer_ids = list(AllChem.EmbedMultipleConfs(
      mol, numConfs=args.num_conformers, params=params))
  if not conformer_ids:
    return None, []

  force_field, results = optimize(mol, args)
  records = []
  for conformer_id, (status, energy) in zip(conformer_ids, results):
    records.append((conformer_id, status, energy))
  records.sort(key=lambda record: (record[2] is None, record[2] or 0.0))

  if args.remove_hydrogens:
    mol = Chem.RemoveHs(mol)
  return mol, (force_field, records)


def main():
  args = parse_args()
  writer = Chem.SDWriter(str(Path(args.output)))
  if writer is None:
    raise SystemExit(f"Cannot open output '{args.output}'")

  molecules_read = 0
  molecules_written = 0
  conformers_written = 0
  try:
    for line_number, smiles, name in read_smiles(args.input):
      molecules_read += 1
      mol = Chem.MolFromSmiles(smiles)
      if mol is None:
        print(f"Invalid SMILES at line {line_number}: {smiles}", file=sys.stderr)
        continue
      mol.SetProp("_Name", name)
      mol.SetProp("INPUT_SMILES", smiles)

      generated, result = conformers(mol, args)
      if generated is None:
        print(f"Embedding failed at line {line_number}: {name}", file=sys.stderr)
        continue

      force_field, records = result
      original_name = name
      for rank, (conformer_id, status, energy) in enumerate(records, 1):
        generated.SetProp("_Name", f"{original_name}_conf_{rank}")
        generated.SetProp("ORIGINAL_NAME", original_name)
        generated.SetIntProp("CONFORMER_RANK", rank)
        generated.SetIntProp("RDKIT_CONFORMER_ID", conformer_id)
        generated.SetProp("EMBEDDING_METHOD", "ETKDGv3")
        generated.SetProp("FORCE_FIELD", force_field)
        optimization_status = {
            0: "converged",
            1: "not_converged",
            2: "force_field_unavailable",
            3: "not_requested",
        }[status]
        generated.SetProp("OPTIMIZATION_STATUS", optimization_status)
        if energy is not None:
          generated.SetDoubleProp("ENERGY_KCAL_MOL", energy)
        elif generated.HasProp("ENERGY_KCAL_MOL"):
          generated.ClearProp("ENERGY_KCAL_MOL")
        writer.write(generated, confId=conformer_id)
        conformers_written += 1
      molecules_written += 1
  finally:
    writer.close()

  print(f"Read {molecules_read} molecules, wrote {conformers_written} conformers "
        f"for {molecules_written} molecules", file=sys.stderr)
  return 0 if conformers_written else 1


if __name__ == "__main__":
  sys.exit(main())

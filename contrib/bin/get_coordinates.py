# Extract coordinates of matched atoms

from absl import app
from absl import flags
from absl import logging

from lillymol import *
from lillymol_io import *
from lillymol_query import *

FLAGS = flags.FLAGS

flags.DEFINE_list('s', [], "Smarts for matched atoms")
flags.DEFINE_list('q', [], "Query files for matched atoms")

def get_coordinate(mol: Molecule, queries: [SubstructureQuery]):
  """Get coordinates of matched atoms in `mol`.
    For each query in `queries` find the matched atoms.
    For each set of matched atoms, write the coordinates of the matched atoms.
    Not sure what the output format should be
    Currently
    name query_number matched_atom_number atom_number asymbol x y z
  """
  got_match = False
  for query_number, query in enumerate(queries):
    matches = query.substructure_search_matches(mol)
    if not matches:
      continue

    got_match = True

    for match in matches:
      for ndx, atom in enumerate(match):
        a = mol[atom]  # the atom object
        print(f"{mol.name()} {query_number} {ndx} {atom} {a.atomic_symbol()} {a.x():.4f} {a.y():.4f} {a.z():.4f}")

  if not got_match:
    return

def main(argv):
  """Temporary tool to fix broken multi-fragment unique smiles problem
  """
  queries = []
  for smarts in FLAGS.s:
    qry = SubstructureQuery()
    if not qry.build_from_smarts(smarts):
      logging.error("Invalid smarts %s", smarts)
      return 1
    queries.append(qry)

  for fname in FLAGS.q:
    qry = SubstructureQuery()
    if not qry.read_proto(fname):
      logging.error("Cannot read proto query %s", fname)
      return 1
    queries.append(qry)

  if not queries:
    logging.error("Must specify one or more queries via the -s or -q options")
    return 1

  if len(argv) != 2:
    logging.error("Must specify input file as a positional argument")
    return 1

  with ReaderContext(argv[1]) as reader:
    for mol in reader:
      get_coordinate(mol, queries)

if __name__ == '__main__':
  app.run(main)

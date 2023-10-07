""" Collate multiple fragstat proto files generated by dicer.
When run with the -B fragstatproto option, dicer generates files
of Dicer.DicerFragment protos.
This script accumulates any number of these and writes a combined
result.
The output file can go directly into

iwbdb_load.sh -v -c 1 -d name.bdb <generated here>
"""

import sys
from typing import Dict

from google.protobuf import text_format

from absl import app
from absl import flags
from absl import logging

import dicer_fragments_pb2

FLAGS = flags.FLAGS

flags.DEFINE_boolean('verbose', False, 'Verbose output')

# Currently has no options

def main(argv):
  """Accumulate the files in ARGV to a single output
  """
  if len(argv) == 1:
    print('Must specify files to be processed', file=sys.stderr)
    exit(1)

  # A mapping from unique smiles to DicerFragment.
  seen: Dict[str, dicer_fragments_pb2.DicerFragment] = {}
  for fname in argv[1:]:
    if FLAGS.verbose:
      logging.info('Processing %s', fname)
    with open(fname, 'r') as reader:
      for line in reader.readlines():
        smiles, rest = line.split(maxsplit=1)
        proto = text_format.Parse(rest, dicer_fragments_pb2.DicerFragment())
        if proto.smi in seen:
          seen[proto.smi].n += 1
        else:
          seen[proto.smi] = proto

  logging.info('Have data on %d fragments', len(seen))
  max_count = 0
  max_atom_count = 0
  for smi, proto in seen.items():
    print(f'{smi} {text_format.MessageToString(proto, as_one_line=True)}')
    if proto.n > max_count:
      max_count = proto.n
    if proto.nat > max_atom_count:
      max_atom_count = proto.nat

  logging.info('Max found %d, max atom count %d', max_count, max_atom_count)
        

if __name__ == '__main__':
  app.run(main)
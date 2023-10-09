"""Process textproto output files from tp1_pipe
"""

from collections import defaultdict # type: ignore
from typing import Dict

from absl import app # type: ignore
from absl import logging # type: ignore

from google.protobuf import text_format

from Molecule_Tools import demerit_pb2

def accumulate_reasons(proto: demerit_pb2.Molecule,
                       reasons: Dict[str, int]):
  """Demo app of what might be done with proto output from the medchem rules.
    
  Args:
    proto: proto containing info on rules that hit a molecule.
    reasons: a dictionary from rule to count. We will update this depending
      on the rules found in `proto`
  """

  for query_match in proto.query_match:
    reasons[query_match.name] += 1

def medchem_rules_summary(argv):
  """Demo app to show reading and parsing of textproto output from tp1_pipe.
  """

  if len(argv) == 1:
    logging.info("Must specify one or more input files\n")
    return 1

  reasons = defaultdict(int)
  molecules_read = 0
  rejected = 0
  no_demerits = 0
  for fname in argv[1:]:
    logging.info("Opening %s", fname)
    with open(fname, 'r') as reader:
      for line in reader:
        molecules_read += 1
        message = text_format.Parse(line, demerit_pb2.Molecule())
        # Here do whatever is to be done with each individual outcome.
        accumulate_reasons(message, reasons)
        if message.rejected:
          rejected += 1
        if len(message.query_match) == 0:
          no_demerits += 1

  for reason,count in reasons.items():
    print(f"{count} instances of {reason}")
  print(f'Read {molecules_read} molecules, {rejected} rejected and {no_demerits} no demerits')

  return 0

if __name__ == '__main__':
  app.run(medchem_rules_summary)

"""Generate random, recursive smarts for benzene.
Generates horrible recursive smarts for testing smarts matchers.
"""
import numpy as np
import re
import sys
from typing import List, Text

from absl import app
from absl import flags

FLAGS = flags.FLAGS
flags.DEFINE_integer("nsmarts", 100, "Number of smarts to generate")
flags.DEFINE_bool("write_last", False, "Only write the last smarts")

# Possible smarts for one of the carbons in benzene.
atoms = [
    'c', '[ch]', '[ac]', '[acr6]', '[cR1]', '[#6aH]', '[c^n]', '[cX3H]',
    '[x2Hc]', '[0#6+0a]', '[cH1a]', '[$(cc)]', 'c', '[0#6H1r6R1D2X3x2+0a]',
    '[cD2]',
]

# We assemble a list of the tokens in the smarts, and randomly select one to replace.
# These tokens, cannot be replaced.
OPEN_GROUP = '[$('
CLOSE_GROUP = ')]'


def recursive_smarts(atoms: List[Text]) -> str:
  """Generate a random smarts for benzene"""
  ring = np.random.randint(1, 99)
  if ring > 9:
    ring = f'%{ring}'
  return [OPEN_GROUP, atoms[0], str(ring), *atoms[1:6], str(ring), CLOSE_GROUP]


def generate_smarts(argv):
  """Generate recursive smarts for benzene"""
  del argv

  # The starting smarts
  xsmarts = recursive_smarts(['c', 'c', 'c', 'c', 'c', 'c'])

  produced = 0
  to_produce = FLAGS.nsmarts
  only_write_last = FLAGS.write_last
  longest = 0
  while produced < to_produce:

    np.random.shuffle(atoms)
    smarts = recursive_smarts(atoms[0:6])

    pos = np.random.randint(0, len(xsmarts))
    if xsmarts[pos] == OPEN_GROUP or xsmarts[pos] == CLOSE_GROUP or re.search(
        r'^%*\d{1,2}$', xsmarts[pos]):
      continue

    del xsmarts[pos]   # Remove old.
    xsmarts[pos:pos] = smarts    # Insert new.

    produced += 1
    if only_write_last and produced < to_produce:
      continue

    smiles = ''.join(xsmarts)
    if len(smiles) > longest:
      longest = len(smiles)

    print(f'{smiles} #{produced}')
    # If using tsubstructure, test with
    # tsubstructure -q S:/path/to/this/output benzene.smi
    # I have observed this working with nsmarts=500, but much after that
    # dynamic memory allocation of the queries gets too large.
    # If only writing the last query, I have observed nsmarts=50000 working,
    # the resulting smarts was 2.0MB.

  print(f'Longest smiles {longest} bytes', file=sys.stderr)

if __name__ == '__main__':
  app.run(generate_smarts)

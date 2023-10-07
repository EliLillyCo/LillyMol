"""Demonstrate python processing of a serialized nnbr::ListOfNearNeighbours proto"""

from statistics import mean

from absl import app
from absl import flags
from absl import logging

from Utilities.GFP_Tools import nearneighbours_pb2

def main(argv:str):
  """Read and process a binary nnbr::ListOfNearNeighbours.
  Args:
    argv: file containing serialized proto.
  """
  if len(argv) != 2:
    logging.fatal('Must specify input file as an argument')

  # Make an empty proto.
  proto = nearneighbours_pb2.ListOfNearNeighbours()
  logging.info('Reading from %s', argv[1])
  # Read and parse the file.
  with open(argv[1], "rb") as reader:
    proto.ParseFromString(reader.read())

  # The ListOfNearNeighbours message has a repeated field nearneighbours, one for
  # each needle. Report the number of needles in the proto.
  logging.info('Read data on %d needles', len(proto.nearneighbours))

  # For each needle, report the id of the needle, the number of neighbours
  # it has and stats on the distances.
  for needle in proto.nearneighbours:
    distances = [nbr.dist for nbr in needle.nbr]
    print(f'{needle.name} has {len(needle.nbr)} nbrs closest {distances[0]:.4f} ave {mean(distances):.4f}')
      

if __name__ == '__main__':
  app.run(main)

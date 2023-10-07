# Read a data file and apply min/max scaling to a
# selected column, writing the min/max values to a proto.
# Ultimately not used because too slow compared to C++, and
# would need to implement unscaling and have the ability to
# read from a pipe.

import csv

from absl import app
from absl import flags
from google.protobuf.json_format import MessageToJson
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

from Utilities.General import feature_scaling_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("input", None, "Input data file to process")
flags.DEFINE_string("output", None, "proto file that will contain the range")
flags.DEFINE_string("sep", ' ', "Input file token delimiter")
flags.DEFINE_integer("header", 1, "Number of header records")
flags.DEFINE_integer("col", 2, "The column to process")
flags.DEFINE_string("proto", None, "Name of FeatureScaling proto to create")

def feature_scaling(unused_argv):
  """Apply min/max scaling to 
  """
  del unused_argv

  col = FLAGS.col - 1

  data = pd.read_csv(FLAGS.input, sep=FLAGS.sep, header=FLAGS.header, usecols=[0, col])

  mycol = np.array(data.iloc[:,col])
  mymin = np.min(mycol)
  mymax = np.max(mycol)
  data.iloc[:,col] = (mycol - mymin) / (mymax - mymin)
  data.to_csv(FLAGS.output, FLAGS.sep, index=False)

  if not FLAGS.proto:
    return

  proto = feature_scaling_pb2.FeatureScaling()
  proto.min = mymin
  proto.max = mymax
  proto.nsamples = mycol.shape[0]
  proto.mean = np.mean(mycol)
  with open(FLAGS.proto, "wb") as output:
    output.write(proto.SerializeToString())
  print(MessageToJson(proto))

if __name__ == '__main__':
  flags.mark_flag_as_required('input')
  flags.mark_flag_as_required('output')
  app.run(feature_scaling)

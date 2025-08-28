# Evaluate an xgboost descriptor model built with xgboost_make

import os
import sys

import joblib
import pandas as pd
import random_forest_model_pb2
from absl import app, flags, logging
from google.protobuf import text_format

from pytk.xgbd.linear_scaling_pb2 import LinearScalingData

FLAGS = flags.FLAGS

flags.DEFINE_string("mdir", "", "Model directory")

def read_rescaling_data(fname: str):
  """Read a linear_scaling_pb2.LinearScalingData serialized proto from `fname`.
  """
  with open(fname, "rb") as input:
    data = input.read()

  result = LinearScalingData()
  result.ParseFromString(data)

  return result

def get_model(mdir: str)->tuple:
  """Look for 'model_metadata.txt` in `mdir` and make sure it is OK
    Return a model instantiated from mdir/random_forest.joblib and the
    name of the response.
    If rescaling.dat is present, also return it.
  """
  failed = None, None, None
  fname = os.path.join(mdir, "model_metadata.txt")
  if not os.path.exists(fname):
    logging.error("%s not found", fname)
    return failed

  with open(fname, "r") as reader:
    text = reader.read()

  proto = text_format.Parse(text, random_forest_model_pb2.RandomForestModel())
  if not proto:
    logging.error("Cannot interpret as proto %s", text)
    return failed

  if not proto.response:
    logging.error("No response in %s", fname)
    return failed

  model_file = os.path.join(mdir, "random_forest.joblib")
  if not os.path.exists(model_file):
    logging.error("%s not found", model_file)
    return failed

  model = joblib.load(model_file)

  rescaling = os.path.join(mdir, 'rescaling.dat')
  if os.path.exists(rescaling):
    True
  else:
    rescaling = None

  return model, proto.response, rescaling

def random_forest_evaluate(mdir: str, fname: str)->bool:
  """Read `fname` as descriptors for a model in `mdir`
  """
  if not os.path.isdir(mdir):
    logging.error("Model directory %s not found", mdir)
    return False

  model, response, rescaling_fname = get_model(mdir)
  if not model:
    logging.error("Invalid mode in %s", mdir)
    return False

  if rescaling_fname is not None:
    rescaling = read_rescaling_data(rescaling_fname)
    # print(rescaling, file=sys.stderr)
  else:
    rescaling = None

  data = pd.read_csv(fname, sep=' ', header=0)

  logging.info("Evaluating %d rows", len(data))
  results = model.predict(data.iloc[:,1:])
  print(f"Id RF_{response}")
  if rescaling is None:
    for i in range(len(results)):
      print(f"{data.iloc[i,0]} {results[i]:.4f}")
  else:
    for i in range(len(results)):
      unscaled = (results[i] - rescaling.intercept) / rescaling.slope
      print(f"{data.iloc[i,0]} {unscaled:.4f}")

  return True

def main(argv):
  """Evaluate a random forest descriptor model.
  """
  if len(argv) == 1:
    logging.error("Must specify descriptor file as argument")
    return 1

  if not FLAGS.mdir:
    logging.error("must specify model directory via the --mdir option")
    return 1

  return random_forest_evaluate(FLAGS.mdir, argv[1])

if __name__ == '__main__':
  app.run(main)

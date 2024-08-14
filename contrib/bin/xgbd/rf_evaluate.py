# Evaluate an xgboost descriptor model built with xgboost_make

import os
import re

import joblib
import pandas as pd

from absl import app
from absl import flags
from absl import logging
from google.protobuf import text_format

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor

import random_forest_model_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("mdir", "", "Model directory")

def get_model(mdir: str)->tuple:
  """Look for 'what_kind_of_model` in `mdir` and make sure it is OK
    Return a model instantiated from mdir/random_forest.joblib and the
    name of the response
  """
  fname = os.path.join(mdir, "model_metadata.txt")
  if not os.path.exists(fname):
    logging.error("%s not found", fname)
    return None, None

  with open(fname, "r") as reader:
    text = reader.read()

  proto = text_format.Parse(text, random_forest_model_pb2.RandomForestModel())
  if not proto:
    logging.error("Cannot interpret as proto %s", text)
    return None, None

  if not proto.response:
    logging.error("No response in %s", fname)
    return None, None

  model_file = os.path.join(mdir, "random_forest.joblib")
  if not os.path.exists(model_file):
    logging.error("%s not found", model_file)
    return None

  model = joblib.load(model_file)

  return model, proto.response

def random_forest_evaluate(mdir: str, fname: str)->bool:
  """Read `fname` as descriptors for a model in `mdir`
  """
  if not os.path.isdir(mdir):
    logging.error("Model directory %s not found", mdir)
    return False

  model, response = get_model(mdir)
  if not model:
    logging.error("Invalid mode in %s", mdir)
    return False

  data = pd.read_csv(fname, sep=' ', header=0)

  logging.info("Evaluating %d rows", len(data))
  results = model.predict(data.iloc[:,1:])
  print(f"Id RF_{response}")
  for i in range(len(results)):
    print(f"{data.iloc[i,0]} {results[i]:.4f}")

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

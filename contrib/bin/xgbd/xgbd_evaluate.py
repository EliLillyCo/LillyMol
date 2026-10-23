# Evaluate an xgboost descriptor model built with xgboost_make

import os
import re

import pandas as pd
import xgboost_model_pb2
from absl import app, flags, logging
from google.protobuf import text_format
from xgboost import XGBClassifier, XGBRegressor
from xgbd import class_label_translation_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("mdir", "", "Model directory")

def get_model(mdir: str)->tuple:
  """Look for 'what_kind_of_model` in `mdir` and make sure it is OK
    Return a model instantiated from mdir/xgboost.json and the
    name of the response
  """
  failed = [None, None, None]

  fname = os.path.join(mdir, "model_metadata.txt")
  if not os.path.exists(fname):
    logging.error("%s not found", fname)
    return failed

  with open(fname, "r") as reader:
    text = reader.read()

  proto = text_format.Parse(text, xgboost_model_pb2.XGBoostModel())
  if not proto:
    logging.error("Cannot interpret as proto %s", text)
    return failed

  if not proto.response:
    logging.error("No response in %s", fname)
    return failed

  model_file = os.path.join(mdir, "xgboost.json")
  if not os.path.exists(model_file):
    logging.error("%s not found", model_file)
    return failed

  if proto.classification:
    model = XGBClassifier()
  else:
    model = XGBRegressor()

  model.load_model(model_file)

  return proto, model, proto.response

def read_class_label_translation(mdir:string)->bool:
  """Read the ClassLabelTransation proto in `mdir`.
  """
  fname = os.path.join(mdir, 'class_label_translation.dat')
  with open(fname, "rb") as input:
    serialised = input.read()

  proto = class_label_translation_pb2.ClassLabelTranslation()
  proto.ParseFromString(serialised)

  return proto

def to_class(labels:List[string], value):
  """Given a list of class labels and a model prediction, return the label for that score.
  """
  if value <= 0.5:
    return labels[0]
  else:
    return labels[1]
 


def xgboost_evaluate(mdir: str, fname: str)->bool:
  """Read `fname` as descriptors for a model in `mdir`
  """
  if not os.path.isdir(mdir):
    logging.error("Model directory %s not found", mdir)
    return False

  proto, model, response = get_model(mdir)
  if not model:
    logging.error("Invalid mode in %s", mdir)
    return False

  print(f"classification {proto.classification}")
  if proto.classification:
    xref = read_class_label_translation(mdir)
    classes = [False] * 2
    print(f"xref {xref}")
    for k, v in enumerate(xref.to_numeric):
      print(f"k #{k} v #{v}")
      classes[k] = v

  data = pd.read_csv(fname, sep=' ', header=0)

  logging.info("Evaluating %d rows", len(data))
  results = model.predict(data.iloc[:,1:])
  print(f"Id XGBD_{response}")
  if proto.classification:
    for i in range(len(results)):
      print(f"{data.iloc[i,0]} {to_class(classes, results[i]):.4f}")
  else:
    for i in range(len(results)):
      print(f"{data.iloc[i,0]} {results[i]:.4f}")

  return True

def main(argv):
  """Evaluate an xgboost descriptor model.
  """
  if len(argv) == 1:
    logging.error("Must specify descriptor file as argument")
    return 1

  if not FLAGS.mdir:
    logging.error("must specify model directory via the --mdir option")
    return 1


  return xgboost_evaluate(FLAGS.mdir, argv[1])

if __name__ == '__main__':
  app.run(main)

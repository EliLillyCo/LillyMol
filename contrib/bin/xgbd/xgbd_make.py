# Build and commit an xgboost model
# Deliberately simplistic in approach

import os

import pandas as pd
import sklearn
from matplotlib import pyplot

from xgboost import plot_importance
from xgboost import XGBClassifier
from xgboost import XGBRegressor

from absl import app
from absl import flags
from absl import logging
from google.protobuf import text_format

import xgboost_model_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("activity", "", "Name of training set activity file")
flags.DEFINE_string("desc", "", "Name of training set descriptor file")
flags.DEFINE_boolean("classification", False, "True if this is a classification task")
flags.DEFINE_string("mdir", "", "Directory into which the model is placed")
flags.DEFINE_integer("max_num_features", 0, "Maximum number of features to plot in variable importance")
flags.DEFINE_string("feature_importance", "", "File containing feature importance values")
flags.DEFINE_integer("xgverbosity", 0, "xgboost verbosity")
flags.DEFINE_string("proto", "", "A file containing an XGBoostParameters proto")
flags.DEFINE_float("eta", 0.3, "xgboost learning rate parameter eta")
flags.DEFINE_integer("max_depth", 6, "xgboost max depth")
flags.DEFINE_integer("n_estimators", 100, "xboost number of estimators")

class Options:
  def __init__(self):
    self.classification = False
    self.mdir: str = ""
    self.max_num_features: int = 10
    self.verbosity = 0
    self.proto = xgboost_model_pb2.XGBoostParameters()

  def read_proto(self, fname)->bool:
    """Read self.proto from `fname`
    """
    with open(fname, "r") as reader:
      text = reader.read()

    self.proto = text_format.Parse(text, xgboost_model_pb2.XGBoostParameters())
    if not self.proto:
      logging.error("Cannot intpret %s", text)
      return False

    return True
def classification(x, y, options: Options)->bool:
  """build a classification model
  """
  booster = XGBClassifier(verbosity=options.verbosity)
  booster.fit(x, y)

def regression(x, y, options: Options):
  """build a regression model.
  """
  booster = XGBRegressor(verbosity=options.verbosity,
                eta=options.proto.eta,
                max_depth=options.proto.max_depth,
                n_estimators = options.proto.n_estimators)
  booster.fit(x, y)

  booster.save_model(os.path.join(options.mdir, "xgboost.json"))
  if options.max_num_features:
    plot_importance(booster, max_num_features=options.max_num_features)
    pyplot.show()
  if options.feature_importance:
    feature_importance = booster.get_booster().get_score(importance_type='weight')
    feature_importance = sorted(feature_importance.items(), key=lambda x:x[1])
    if options.feature_importance:
      with open(os.path.join(options.mdir, options.feature_importance), "w") as writer:
        # Write a markdown table, easy to undo if needed.
        print("| Feature | Weight |", file=writer)
        print("| ------- | ------ |", file=writer)
        for f, i in feature_importance:
          print(f"| {f} | {i} |", file=writer)

  # config = booster.save_config()

  return True

def build_xgboost_model(descriptor_fname: str,
                        activity_fname: str,
                        options: Options)->bool:
  """Build an xgboost model on the data in `descriptor_fname` and
     `activity_fname`.
    This function does data preprocessing.
  """

  descriptors = pd.read_csv(descriptor_fname, sep=' ', header=0, low_memory=False)
  logging.info("Read %d rows and %d columns from %s", len(descriptors),
                descriptors.shape[1], descriptor_fname)
  activity = pd.read_csv(activity_fname, sep=' ', header=0)
  logging.info("Read %d rows from %s", activity.shape[0], activity_fname)


  descriptors.rename(columns={descriptors.columns[0]: "Id"}, inplace=True)
  activity.rename(columns={activity.columns[0]: "Id"}, inplace=True)
  combined = pd.concat([activity.set_index("Id"),
                        descriptors.set_index("Id")], axis=1, join='inner').reset_index() 
  if len(combined) != len(descriptors):
    logging.error("Combined set has %d rows", len(combined))
    return 1

  if not os.path.isdir(options.mdir):
    os.mkdir(options.mdir)

  y = combined.iloc[:,1].to_numpy()

  x = combined.iloc[:,2:]
  x.apply(pd.to_numeric).to_numpy()

  rc = False
  if options.classification:
    rc = classification(x, y, options)
  else:
    rc = regression(x, y, options)

  if not rc:
    return False

  response = activity.columns[1]

  proto = xgboost_model_pb2.XGBoostModel();
  proto.model_type = "XGBD"
  proto.classification = False
  proto.response = response
  proto.parameters.CopyFrom(options.proto)
  with open(os.path.join(options.mdir, "model_metadata.txt"), "w") as f:
    f.write(text_format.MessageToString(proto))

  return True

def main(argv):
  """Build xgboost models from activity file and descriptor file.
  """
  if not FLAGS.activity:
    logging.error("Must specifythe name of the activity file with the --activity option")
    return False
  if not FLAGS.desc:
    logging.error("Must specifythe name of the descriptor file with the --desc option")
    return False
  if not FLAGS.mdir:
    logging.error("Must specifyi the model directory via the --mdir option")
    return False

  options = Options()
  options.classification = FLAGS.classification
  options.mdir = FLAGS.mdir
  options.max_num_features = FLAGS.max_num_features
  options.feature_importance = FLAGS.feature_importance
  options.verbosity = FLAGS.xgverbosity

  # Build the proto first, and then anything that might overwrite it.
  if FLAGS.proto:
    if not options.read_proto(FLAGS.proto):
      logging.error("Cannot read textproto parameters %s", FLAGS.proto)
      return False
  else:
    options.proto.eta = FLAGS.eta
    options.proto.max_depth = FLAGS.max_depth
    options.proto.n_estimators = FLAGS.n_estimators

  if not build_xgboost_model(FLAGS.desc, FLAGS.activity, options):
    logging.error("Model %s not build", options.mdir)
    return False

  return True

if __name__ == '__main__':
  app.run(main)

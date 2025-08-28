# Build and commit an xgboost model
# Deliberately simplistic in approach

import os
import subprocess
import sys
from typing import List

import pandas as pd
import xgboost_model_pb2
from absl import app, flags, logging
from google.protobuf import text_format
from matplotlib import pyplot
from xgboost import XGBClassifier, XGBRegressor, plot_importance

FLAGS = flags.FLAGS

flags.DEFINE_string("activity", "", "Name of training set activity file")
flags.DEFINE_boolean("classification", False, "True if this is a classification task")
flags.DEFINE_string("mdir", "", "Directory into which the model is placed")
flags.DEFINE_integer("max_num_features", 0, "Maximum number of features to plot in variable importance")
flags.DEFINE_string("feature_importance", "", "Compute feature importance. Use 'def' to use default file")
flags.DEFINE_integer("xgverbosity", 0, "xgboost verbosity")
flags.DEFINE_string("proto", "", "A file containing an XGBoostParameters proto")
flags.DEFINE_float("eta", 0.4, "xgboost learning rate parameter eta")
flags.DEFINE_integer("max_depth", 5, "xgboost max depth")
flags.DEFINE_integer("n_estimators", 500, "xboost number of estimators")
flags.DEFINE_float("subsample", 1.0, "subsample ratio for training instances")
flags.DEFINE_float("colsample_bytree", 1.0, "subsampling occurs once for every tree constructed")
flags.DEFINE_float("colsample_bylevel", 1.0, "subsampling occurs once for every new depth level reached")
flags.DEFINE_float("colsample_bynode", 1.0, "subsampling occurs once for every time a new split is evaluated")
flags.DEFINE_enum("tree_method", "auto", ["auto", "exact", "approx", "hist"], "tree construction method: auto exact approx hist")
flags.DEFINE_integer("nthreads", 8, "number of threads to use, default is 8")
flags.DEFINE_boolean("rescore", False, "Rescore the training set to establish linear correction function")


class Options:
  def __init__(self):
    self.classification = False
    self.mdir: str = ""
    self.max_num_features: int = 10
    self.descriptor_fname: str = ""
    self.activity_fname: str = ""
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

def to_array(input:str) -> List[str]:
  return input.split(' ')

def classification(x, y, options: Options)->bool:
  """build a classification model
    Args:
      x: feature matrix
      y: response - must be translated to 0,1. Not implemented...
  """
  booster = XGBClassifier(verbosity=options.verbosity)
  booster.fit(x, y)

  booster.save_model(os.path.join(options.mdir, "xgboost.json"))

def regression(x, y, options: Options):
  """build a regression model.
  """
  match options.proto.tree_method:
    case xgboost_model_pb2.AUTO:
      tree_method = 'auto'
    case xgboost_model_pb2.EXACT:
      tree_method = 'exact'
    case xgboost_model_pb2.APPROX:
      tree_method = 'approx'
    case xgboost_model_pb2.HIST:
      tree_method = 'hist'
    case _:
      tree_method = 'hist'

  nthread = None
  if FLAGS.nthreads > 0:
    nthread = FLAGS.nthreads

  booster = XGBRegressor(verbosity=options.verbosity,
                eta=options.proto.eta,
                max_depth=options.proto.max_depth,
                n_estimators=options.proto.n_estimators,
                colsample_bytree=options.proto.colsample_bytree,
                colsample_bylevel=options.proto.colsample_bylevel,
                colsample_bynode=options.proto.colsample_bynode,
                subsample=options.proto.subsample,
                tree_method=tree_method,
                nthread=nthread
                )

  booster.fit(x, y)

  booster.save_model(os.path.join(options.mdir, "xgboost.json"))
  logging.info("Saved model to %s", os.path.join(options.mdir, "xgboost.json"))

  if options.max_num_features:
    plot_importance(booster, max_num_features=options.max_num_features)
    pyplot.show()

  if len(options.feature_importance) > 0:
    for itype in ["weight", "gain", "cover"]:
      feature_importance = booster.get_booster().get_score(importance_type=itype)
      feature_importance = sorted(feature_importance.items(), key=lambda x:x[1], reverse=True)

      fname = options.feature_importance
      if fname == "def" or fname == "DEF":
        fname = os.path.join(options.mdir, f'feature_importance.{itype}.txt')
      else:
        # Or should this be placed in `mdir` by default?
        fname = f"{fname}.{itype}.txt"

      with open(fname, "w") as writer:
        print("Feature Weight", file=writer)
        for f, i in feature_importance:
          print(f"{f} {i}", file=writer)

  return True


def build_xgboost_model(descriptor_fname: str,
                        activity_fname: str,
                        options: Options)->bool:
  """Build an xgboost model on the data in `descriptor_fname` and
     `activity_fname`.
    This function does data preprocessing.
  """

  descriptors = pd.read_csv(descriptor_fname, sep=' ', header=0, low_memory=False, na_values=['.'])
  logging.info("Read %d rows and %d columns from %s", len(descriptors),
                descriptors.shape[1], descriptor_fname)
  activity = pd.read_csv(activity_fname, sep=' ', header=0)
  logging.info("Read %d rows from %s", activity.shape[0], activity_fname)


  descriptors.rename(columns={descriptors.columns[0]: "Name"}, inplace=True)
  activity.rename(columns={activity.columns[0]: "Name"}, inplace=True)
  combined = pd.concat([activity.set_index("Name"),
                        descriptors.set_index("Name")], axis=1, join='inner').reset_index() 
  if len(combined) != len(descriptors):
    logging.error("Combined set has %d rows, need %d", len(combined), len(descriptors))
    return FALSE

  if not os.path.isdir(options.mdir):
    os.mkdir(options.mdir)

  combined.to_csv(os.path.join(options.mdir, "train.xy"), sep= ' ', index=False)

  y = combined.iloc[:,1].to_numpy()

  x = combined.iloc[:,2:]
  features = x.columns
  x.apply(pd.to_numeric).to_numpy()

  options.descriptor_fname = descriptor_fname
  options.activity_fname = activity_fname

  rc = False
  if options.classification:
    rc = classification(x, y, options)
  else:
    rc = regression(x, y, options)

  if not rc:
    logging.info("Model did not build")
    return False

  response = activity.columns[1]

  proto = xgboost_model_pb2.XGBoostModel()
  proto.model_type = "XGBD"
  proto.classification = options.classification
  proto.response = response
  proto.parameters.CopyFrom(options.proto)

  for (column, feature) in enumerate(features):
    proto.name_to_col[feature] = column

  with open(os.path.join(options.mdir, "model_metadata.txt"), "w") as f:
    f.write(text_format.MessageToString(proto))
  with open(os.path.join(options.mdir, "model_metadata.dat"), "wb") as f:
    f.write(proto.SerializeToString())

  return True

def option_present(flag)->bool:
  """Return true if the option `flag` is in sys.argv
    Args:
      argv: usually the command line
      flag: a command line option. We look for -flag and --flag in argv.
  """
  if '-'+flag in sys.argv:
    return True
  if '--'+flag in sys.argv:
    return True
  return False


def rescore_training_set(options)->bool:
  """ A model has just been build in `options.mdir`.
      Rescore the training set and store the results.
  """
  train_pred = os.path.join(options.mdir, 'train.pred')
  with open(train_pred, 'w') as output:
    cmd = f"xgbd_evaluate.sh -mdir {options.mdir} {options.descriptor_fname}"
    subprocess.run(to_array(cmd), stdout=output, text=True)

  if not os.path.exists(train_pred):
    logging.error("%s did not create %s", cmd, train_pred)
    return False

  train_stats = os.path.join(options.mdir, 'train.stats')
  rescaling = os.path.join(options.mdir, 'rescaling.textproto')
  cmd = f"iwstats.sh -Y allequals -w -E {options.activity_fname} -p 2 -C {rescaling} {train_pred}"
  with open(train_stats, 'w') as output:
    subprocess.run(to_array(cmd), stdout=output, text=True)

  if not os.path.exists(rescaling):
    logging.error("%s did not create %s", cmd, rescaling)
    return False

  return True

def main(argv):
  """Build xgboost models from activity file and descriptor file.
  """
  if not FLAGS.activity:
    logging.error("Must specifythe name of the activity file with the --activity option")
    return False
  if len(argv) == 1:
    logging.error("Must specifythe name of the descriptor file as argument")
    return False
  if not FLAGS.mdir:
    logging.error("Must specify the model directory via the --mdir option")
    return False

  options = Options()
  options.classification = FLAGS.classification
  options.mdir = FLAGS.mdir
  options.max_num_features = FLAGS.max_num_features
  options.feature_importance = FLAGS.feature_importance
  options.verbosity = FLAGS.xgverbosity

  # Build the proto first.
  # After that is done, we check for command line arguments that would
  # over-ride what has come in from the proto.
  if FLAGS.proto:
    if not options.read_proto(FLAGS.proto):
      logging.error("Cannot read textproto parameters %s", FLAGS.proto)
      return False

  # Overrides from the command line.
  # If the proto does not have a value, just use the default from FLAGS.
  if not options.proto.HasField("eta"):
    options.proto.eta = FLAGS.eta

  if not options.proto.HasField("max_depth"):
    options.proto.max_depth = FLAGS.max_depth

  if not options.proto.HasField("n_estimators"):
    options.proto.n_estimators = FLAGS.n_estimators

  if not options.proto.HasField("subsample"):
    options.proto.subsample = FLAGS.subsample

  if not options.proto.HasField("colsample_bytree"):
    options.proto.colsample_bytree = FLAGS.colsample_bytree

  if not options.proto.HasField("colsample_bylevel"):
    options.proto.colsample_bylevel = FLAGS.colsample_bylevel

  if not options.proto.HasField("colsample_bynode"):
    options.proto.colsample_bynode = FLAGS.colsample_bynode

  if option_present("tree_method"):
    match FLAGS.tree_method:
      case "auto":
        options.proto.tree_method = xgboost_model_pb2.AUTO
      case "exact":
        options.proto.tree_method = xgboost_model_pb2.EXACT
      case "approx":
        options.proto.tree_method = xgboost_model_pb2.APPROX
      case "hist":
        options.proto.tree_method = xgboost_model_pb2.HIST
      case _:   # Cannot happen because this is DEFINE_enum
        print(f"Unrecognised tree method {FLAGS.tree_method}", file=sys.stderr)
        return False
  elif not options.proto.HasField("tree_method"):
    options.proto.tree_method = xgboost_model_pb2.AUTO

  if not build_xgboost_model(argv[1], FLAGS.activity, options):
    logging.error("Model %s not build", options.mdir)
    return False

  if FLAGS.rescore:
    rescore_training_set(options)

  # zero return code for success.
  return 0

if __name__ == '__main__':
  app.run(main)

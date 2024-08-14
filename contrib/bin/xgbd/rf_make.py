# Build and commit an RF model
# Deliberately simplistic in approach

import os

import joblib
import pandas as pd
import sklearn

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor

from absl import app
from absl import flags
from absl import logging
from google.protobuf import text_format

import random_forest_model_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("activity", "", "Name of training set activity file")
flags.DEFINE_boolean("classification", False, "True if this is a classification task")
flags.DEFINE_string("mdir", "", "Directory into which the model is placed")
flags.DEFINE_integer("max_num_features", 0, "Maximum number of features to plot in variable importance")
flags.DEFINE_boolean("feature_importance", False, "Create mdir/feature_importance.txt or not")
flags.DEFINE_integer("rfverbosity", 0, "RF verbosity")
flags.DEFINE_string("proto", "", "A file containing an RandomForestModel proto")
flags.DEFINE_integer("max_depth", 12, "max depth")
flags.DEFINE_integer("n_estimators", 100, "number of estimators")
flags.DEFINE_integer("min_samples_split", 4, "minum number of samples required to split an internal node")
flags.DEFINE_integer("min_samples_leaf", 2, "minum number of samples required to be in a leaf node")
flags.DEFINE_integer("n_jobs", -1, "parallelism, -1 means use all available processors")

class Options:
  def __init__(self):
    self.classification = False
    self.mdir: str = ""
    self.max_num_features: int = 10
    self.feature_importance: bool = False
    self.proto = random_forest_model_pb2.RandomForestParameters()

  def read_proto(self, fname)->bool:
    """Read self.proto from `fname`
    """
    with open(fname, "r") as reader:
      text = reader.read()

    self.proto = text_format.Parse(text, random_forest_model_pb2.RandomForestParameters())
    if not self.proto:
      logging.error("Cannot intpret %s", text)
      return False

    return True

# Note that we are doing just impurity based feature importance which is known
# to have a bias against features with just a few values. So integer values
# will likely show up as less important than floating point descriptors.
def do_feature_importance(booster, feature_names, options: Options):
  logging.info("Doing feature importance")

  importances = booster.feature_importances_

  forest_importances = pd.Series(importances, index=feature_names).sort_values(ascending=False)

# plotting never worked as I wanted, tabular output seems good enough.
# If ever anyone wants to implement this, need to sort at the end not above...
# std = np.std([tree.feature_importances_ for tree in booster.estimators_], axis=0)
# fig, ax = plt.subplots()
# forest_importances.plot.bar(yerr=std, ax=ax)
# ax.set_title("Feature importances using MDI")
# ax.set_ylabel("Mean decrease in impurity")
# fig.tight_layout()
# fig.savefig(os.path.join(options.mdir, "feature_importance.png"))

  forest_importances = pd.DataFrame(forest_importances).reset_index()
  forest_importances.columns = ["Feature", "Importance"]
  forest_importances.to_csv(os.path.join(options.mdir, "feature_importance.txt"), sep=' ', index=False)

def keywords_from(proto):
  """Transfer keyword arguments from `proto` to a dictionary.
  """
  result = {}
  if proto.n_estimators:
    result['n_estimators'] = proto.n_estimators
  if proto.max_depth:
    result['max_depth'] = proto.max_depth
  if proto.min_samples_leaf:
    result['min_samples_leaf'] = proto.min_samples_leaf
  if proto.min_samples_split:
    result['min_samples_split'] = proto.min_samples_split
  if proto.n_jobs > 0:
    result['n_jobs'] = proto.n_jobs
  if proto.verbose:
    result['verbose'] = True

  return result

def classification(x, y, options: Options):
  """build a classification model
  """
  args = keywords_from(options.proto)

  logging.info("Build model with %s", args)
  booster = RandomForestClassifier(**args)
  booster.fit(x, y)

  logging.info("Saving model to %s", os.path.join(options.mdir, "random_forest.joblib"))
  joblib.dump(booster, os.path.join(options.mdir, "random_forest.joblib"))

  if len(options.feature_importance):
    do_feature_importance(booster, x.columns, options)

  return True

def regression(x, y, options: Options):
  """build a regression model.
  """
  args = keywords_from(options.proto)
  
  logging.info("Build model with %s", args)
  booster = RandomForestRegressor(**args)
  booster.fit(x, y)

  logging.info("Saving model to %s", os.path.join(options.mdir, "random_forest.joblib"))
  joblib.dump(booster, os.path.join(options.mdir, "random_forest.joblib"))

  if options.feature_importance:
    do_feature_importance(booster, x.columns, options)

  return True

def build_random_forest_model(descriptor_fname: str,
                        activity_fname: str,
                        options: Options)->bool:
  """Build a random forest model on the data in `descriptor_fname` and
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
    logging.error("Combined set has %d rows, need %d", len(combined), len(descriptors))
    return 1

  if not os.path.isdir(options.mdir):
    os.mkdir(options.mdir)

  y = combined.iloc[:,1].to_numpy()

  x = combined.iloc[:,2:]
  features = x.columns
  x.apply(pd.to_numeric).to_numpy()

  rc = False
  if options.classification:
    rc = classification(x, y, options)
  else:
    rc = regression(x, y, options)

  if not rc:
    return False

  response = activity.columns[1]

  proto = random_forest_model_pb2.RandomForestModel();
  proto.model_type = "RF"
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

def main(argv):
  """Build random_forest models from activity file and descriptor file.
  """
  if not FLAGS.activity:
    logging.error("Must specifythe name of the activity file with the --activity option")
    return False
  if not FLAGS.mdir:
    logging.error("Must specify the model directory via the --mdir option")
    return False
  # Descriptor file must be the only command line argument
  if len(argv) == 1:
    logging.error("Must specify the descriptor as a command line argument")
    return False

  options = Options()
  options.classification = FLAGS.classification
  options.feature_importance = FLAGS.feature_importance
  options.mdir = FLAGS.mdir
  options.verbosity = FLAGS.rfverbosity

  # Build the proto first, and then anything that might overwrite it.
  if FLAGS.proto:
    if not options.read_proto(FLAGS.proto):
      logging.error("Cannot read textproto parameters %s", FLAGS.proto)
      return 1

  if FLAGS.n_estimators:
    options.proto.n_estimators = FLAGS.n_estimators
  if FLAGS.max_depth:
    options.proto.max_depth = FLAGS.max_depth
  if FLAGS.min_samples_split:
    options.proto.min_samples_split = FLAGS.min_samples_split
  if FLAGS.min_samples_leaf:
    options.proto.min_samples_leaf = FLAGS.min_samples_leaf
  if FLAGS.n_jobs > 0:
    options.proto.n_jobs = FLAGS.n_jobs
  options.proto.verbose = FLAGS.rfverbosity

  # print(options.proto)

  if not build_random_forest_model(argv[1], FLAGS.activity, options):
    logging.error("Model %s not build", options.mdir)
    return 1

  # zero return code for success.
  return 0

if __name__ == '__main__':
  app.run(main)

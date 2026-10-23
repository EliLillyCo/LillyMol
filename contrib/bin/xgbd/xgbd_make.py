# Build and commit an xgboost model
# Deliberately simplistic in approach

import os
import subprocess

import numpy as np
import pandas as pd

from absl import app, flags, logging
from google.protobuf import text_format
from google.protobuf import json_format
from matplotlib import pyplot
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from xgboost import XGBClassifier, XGBRegressor, plot_importance
from xgboost.core import XGBoostError
from xgbd import xgboost_model_pb2
from xgbd import class_label_translation_pb2

FLAGS = flags.FLAGS

flags.DEFINE_integer("min_points", 10, "do NOT build a model if there are fewer than min_points")
flags.DEFINE_string("activity", "", "Name of training set activity file")
flags.DEFINE_boolean("classification", False, "True if this is a classification task")
flags.DEFINE_string("mdir", "", "Directory into which the model is placed")
flags.DEFINE_integer("max_num_features", 0, "Maximum number of features to plot in variable importance")
flags.DEFINE_string("feature_importance", "", "Compute feature importance. Use 'def' to use default file")
flags.DEFINE_boolean("optuna", False, "Optimize regression parameters with Optuna")
flags.DEFINE_integer("xgverbosity", 0, "xgboost verbosity")
flags.DEFINE_string("proto", "", "A file containing an XGBoostParameters proto")
flags.DEFINE_float("eta", 0.4, "xgboost learning rate parameter eta")
flags.DEFINE_integer("max_depth", 5, "xgboost max depth")
flags.DEFINE_integer("n_estimators", 500, "xboost number of estimators")
flags.DEFINE_float("subsample", 1.0, "subsample ratio for training instances")
flags.DEFINE_float("min_child_weight", 1.0, "Minimum sum of instance weights in a child")
flags.DEFINE_float("colsample_bytree", 1.0, "subsampling occurs once for every tree constructed")
flags.DEFINE_float("colsample_bylevel", 1.0, "subsampling occurs once for every new depth level reached")
flags.DEFINE_float("colsample_bynode", 1.0, "subsampling occurs once for every time a new split is evaluated")
flags.DEFINE_float("reg_lambda", 1.0, "L2 regularization")
flags.DEFINE_float("reg_alpha", 0.0, "L1 regularization")
flags.DEFINE_float("gamma", 0.0, "minimum loss reduction")
flags.DEFINE_enum("tree_method", "auto", ["auto", "exact", "approx", "hist"], "tree construction method: auto exact approx hist")
flags.DEFINE_integer("nthreads", 8, "number of threads to use, default is 8")
flags.DEFINE_boolean("rescore", False, "Rescore the training set to establish linear correction function")


class Options:
  def __init__(self):
    self.min_points = 10
    self.classification = False
    self.mdir: str = ""
    self.max_num_features: int = 0
    self.descriptor_fname: str = ""
    self.activity_fname: str = ""
    self.verbosity = 0
    self.proto = xgboost_model_pb2.XGBoostParameters()
    self.optuna = False
    self.feature_importance = ""
    self.nthreads = 8

  def read_proto(self, fname)->bool:
    """Read self.proto from `fname`
    """
    with open(fname, "r") as reader:
      text = reader.read()

    self.proto = text_format.Parse(text, xgboost_model_pb2.XGBoostParameters())
    return True

def write_class_label_translation(options: Options, categories:np.array, class_counts:np.array)->bool:
  """Write the class label proto to the model directory.
  """
  proto = class_label_translation_pb2.ClassLabelTranslation()
  for number, (label, count) in enumerate(zip(categories, class_counts)):
    proto.to_numeric[str(label)] = number
    proto.class_count[str(label)] = int(count)


  fname = os.path.join(options.mdir, "class_label_translation.dat")
  with open(fname, "wb") as output:
    serialised = proto.SerializeToString()
    output.write(serialised)

  fname = os.path.join(options.mdir, "class_label_translation.json")
  with open(fname, "w") as output:
    output.write(json_format.MessageToJson(proto))

  return True

def classification(x, y, options: Options)->bool:
  """build a classification model
    Args:
      x: feature matrix
      y: response labels; the minority class is encoded as 1.
  """
  categories, counts = np.unique(y, return_counts=True)
  if len(categories) != 2:
    logging.error("Expected two classes, found %d", len(categories))
    return False
  # Keep the minority class positive, and write exactly the mapping used to fit.
  if counts[0] <= counts[1]:
    categories = categories[::-1]
    counts = counts[::-1]
  y = (y == categories[1]).astype(int)
  booster = XGBClassifier(**model_parameters(options))
  booster.fit(x, y)

  booster.save_model(os.path.join(options.mdir, "xgboost.json"))
  write_class_label_translation(options, categories, counts)

  return True

TREE_METHODS = {
    xgboost_model_pb2.AUTO: "auto",
    xgboost_model_pb2.EXACT: "exact",
    xgboost_model_pb2.APPROX: "approx",
    xgboost_model_pb2.HIST: "hist",
}
PARAMETER_NAMES = (
    "eta", "max_depth", "n_estimators", "min_child_weight", "subsample",
    "colsample_bytree", "colsample_bylevel", "colsample_bynode",
    "reg_alpha", "reg_lambda", "gamma",
)


def model_parameters(options):
  """Use the same parameters for classification, regression, and metadata."""
  params = {name: getattr(options.proto, name) for name in PARAMETER_NAMES}
  params.update(verbosity=options.verbosity,
                tree_method=TREE_METHODS[options.proto.tree_method],
                n_jobs=options.nthreads if options.nthreads > 0 else None)
  return params


def apply_parameter_flags(options):
  """Explicit flags override the proto; defaults only fill missing fields."""
  for canonical, alias in (("eta", "learning_rate"), ("gamma", "min_split_loss"),
                           ("reg_alpha", "alpha"), ("reg_lambda", "lambda")):
    if options.proto.HasField(alias):
      setattr(options.proto, canonical, getattr(options.proto, alias))
  for name in PARAMETER_NAMES:
    if FLAGS[name].present or not options.proto.HasField(name):
      setattr(options.proto, name, getattr(FLAGS, name))
  if FLAGS["tree_method"].present or not options.proto.HasField("tree_method"):
    options.proto.tree_method = next(
        key for key, value in TREE_METHODS.items() if value == FLAGS.tree_method)


def regression_with_optuna(x, y, options: Options):
  """Tune on fixed folds, then choose a tree count for fitting all training rows."""
  import optuna
  from optuna.samplers import TPESampler

  values = np.asarray(x)
  y = np.asarray(y)
  random_state = 42
  # Every trial sees identical validation rows, so scores are comparable.
  folds = list(KFold(n_splits=5, shuffle=True,
                    random_state=random_state).split(values))
  base_params = model_parameters(options)
  base_params.pop("eta")  # Use the learning_rate alias during tuning.

  def objective(trial: optuna.Trial) -> float:
      params = {
          **base_params,
          # Core
          "n_estimators": 10_000,  # large; early stopping finds the effective number
          "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
          "max_depth": trial.suggest_int("max_depth", 2, 12),
          "min_child_weight": trial.suggest_float("min_child_weight", 1e-3, 50.0, log=True),
          "subsample": trial.suggest_float("subsample", 0.5, 1.0),
          "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
          "gamma": trial.suggest_float("gamma", 0.0, 10.0),
          "reg_alpha": trial.suggest_float("reg_alpha", 1e-10, 10.0, log=True),
          "reg_lambda": trial.suggest_float("reg_lambda", 1e-10, 100.0, log=True),

          # Tree method (change if you want GPU)
          "tree_method": base_params["tree_method"],

          # Objective / metric
          "objective": "reg:squarederror",
          "eval_metric": "rmse",

          # Repro / speed
          "random_state": random_state,
          "n_jobs": base_params["n_jobs"],
          "early_stopping_rounds": 200,
      }

      fold_rmses = []
      tree_counts = []
      for fold_idx, (tr_idx, va_idx) in enumerate(folds, start=1):
          x_tr, x_va = values[tr_idx], values[va_idx]
          y_tr, y_va = y[tr_idx], y[va_idx]

          model = XGBRegressor(**params)

          # Early stopping: use a validation set within each CV fold
          model.fit(
              x_tr,
              y_tr,
              eval_set=[(x_va, y_va)],
              verbose=False,
          )

          preds = model.predict(x_va)
          rmse = np.sqrt(mean_squared_error(y_va, preds))
          tree_counts.append(model.best_iteration + 1)
          fold_rmses.append(rmse)

          # Let Optuna prune unpromising trials
          trial.report(float(np.mean(fold_rmses)), step=fold_idx)
          if trial.should_prune():
              raise optuna.TrialPruned()

      # Early stopping cannot use the full training set as its validation set.
      # Refit with the mean selected tree count instead.
      trial.set_user_attr("n_estimators", max(1, int(round(np.mean(tree_counts)))))
      return float(np.mean(fold_rmses))

  study = optuna.create_study(
      direction="minimize",
      sampler=TPESampler(seed=random_state),
      pruner=optuna.pruners.MedianPruner(n_startup_trials=10, n_warmup_steps=1),
  )

  n_trials = 180
  timeout = 6000
  study.optimize(objective, n_trials=n_trials, timeout=timeout, show_progress_bar=False)

  logging.info("Best Optuna parameters: %s", study.best_params)
  for name, value in study.best_params.items():
    setattr(options.proto, "eta" if name == "learning_rate" else name, value)
  options.proto.n_estimators = study.best_trial.user_attrs["n_estimators"]

 
def regression(x, y, options: Options):
  """build a regression model.
  """
  if options.optuna:
    regression_with_optuna(x, y, options)
  booster = XGBRegressor(**model_parameters(options))

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

  def read_table(fname):
    # IDs are strings: preserve leading zeros and literal identifiers such as NA.
    first_column = pd.read_csv(fname, sep=r"\s+", nrows=0).columns[0]
    return pd.read_csv(fname, sep=r"\s+", low_memory=False,
                       dtype={first_column: str}, keep_default_na=False)

  descriptors = read_table(descriptor_fname)
  activity = read_table(activity_fname)
  if descriptors.shape[1] < 2 or activity.shape[1] != 2:
    logging.error("Need descriptors and exactly one activity column")
    return False
  descriptors.rename(columns={descriptors.columns[0]: "Name"}, inplace=True)
  activity.rename(columns={activity.columns[0]: "Name"}, inplace=True)
  if activity.columns[1] in descriptors.columns:
    logging.error("Activity column also occurs in descriptors")
    return False
  if descriptors["Name"].isna().any() or activity["Name"].isna().any():
    logging.error("Missing molecule identifiers")
    return False
  # Merge columns directly to avoid fragmenting a wide mixed-type dataframe.
  combined = pd.merge(activity, descriptors, how="inner", on="Name",
                      validate="one_to_one")
  if len(combined) != len(descriptors):
    logging.error("Combined set has %d rows, need %d", len(combined), len(descriptors))
    return False
  if len(combined) < max(options.min_points, 5 if options.optuna else 1):
    logging.error("Not enough rows in training set: %d", len(combined))
    return False

  y = combined.iloc[:, 1]
  if y.isna().any() or y.isin(["", "."]).any():
    logging.error("Missing activity values")
    return False
  if not options.classification:
    y = pd.to_numeric(y)
    if not np.isfinite(y).all():
      logging.error("Non-finite activity values")
      return False
  y = y.to_numpy()
  x = combined.iloc[:, 2:].replace({".": np.nan, "": np.nan}).apply(pd.to_numeric)
  if np.isinf(x.to_numpy()).any():
    logging.error("Infinite descriptor values")
    return False
  features = x.columns
  os.makedirs(options.mdir, exist_ok=True)
  combined.to_csv(os.path.join(options.mdir, "train.xy"), sep=" ", index=False)

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

def rescore_training_set(options)->bool:
  """ A model has just been built in `options.mdir`.
      Rescore the training set and store the results.
  """
  train_pred = os.path.join(options.mdir, 'train.pred')
  with open(train_pred, 'w') as output:
    cmd = ["xgbd_evaluate.sh", "-mdir", options.mdir, options.descriptor_fname]
    subprocess.run(cmd, stdout=output, text=True, check=True)

  if os.path.getsize(train_pred) == 0:
    logging.error("%s did not create %s", cmd, train_pred)
    return False

  train_stats = os.path.join(options.mdir, 'train.stats')
  rescaling = os.path.join(options.mdir, 'rescaling.textproto')
  cmd = ["iwstats.sh", "-Y", "allequals", "-w", "-E", options.activity_fname,
         "-p", "2", "-C", rescaling, train_pred]
  with open(train_stats, 'w') as output:
    subprocess.run(cmd, stdout=output, text=True, check=True)

  if not os.path.exists(rescaling):
    logging.error("%s did not create %s", cmd, rescaling)
    return False

  return True

def build_from_flags(argv):
  """Build xgboost models from activity file and descriptor file.
  """
  if not FLAGS.activity:
    logging.error("Must specify the name of the activity file with the --activity option")
    return False
  if len(argv) != 2:
    logging.error("Must specify the name of the descriptor file as argument")
    return False
  if not FLAGS.mdir:
    logging.error("Must specify the model directory via the --mdir option")
    return False

  options = Options()
  options.min_points = FLAGS.min_points
  options.classification = FLAGS.classification
  options.mdir = FLAGS.mdir
  options.max_num_features = FLAGS.max_num_features
  options.feature_importance = FLAGS.feature_importance
  options.optuna = FLAGS.optuna
  options.verbosity = FLAGS.xgverbosity
  options.nthreads = FLAGS.nthreads

  if FLAGS.classification and (FLAGS.rescore or FLAGS.optuna):
    logging.error("Classification does not support --rescore or --optuna")
    return False

  # Build the proto first.
  # After that is done, we check for command line arguments that would
  # over-ride what has come in from the proto.
  if FLAGS.proto:
    if not options.read_proto(FLAGS.proto):
      logging.error("Cannot read textproto parameters %s", FLAGS.proto)
      return False

  apply_parameter_flags(options)

  if not build_xgboost_model(argv[1], FLAGS.activity, options):
    logging.error("Model %s not built", options.mdir)
    return False

  if FLAGS.rescore and not rescore_training_set(options):
    return False

  return True


def main(argv):
  """Translate operational failures into a nonzero process exit status."""
  try:
    return 0 if build_from_flags(argv) else 1
  except (OSError, ValueError, text_format.ParseError,
          subprocess.CalledProcessError, XGBoostError) as error:
    logging.error("Model build failed: %s", error)
    return 1

if __name__ == '__main__':
  app.run(main)

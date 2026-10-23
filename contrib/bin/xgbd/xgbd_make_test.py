"""Focused regression tests: python3 -m unittest xgbd.xgbd_make_test."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import warnings

import numpy as np
import pandas as pd
from absl.testing import flagsaver
from xgbd import xgbd_make as make
from xgbd import class_label_translation_pb2, xgboost_model_pb2

if not make.FLAGS.is_parsed():
  make.FLAGS(["test"])


class ModelTest(unittest.TestCase):
  def setUp(self):
    self.tmp = tempfile.TemporaryDirectory()
    self.addCleanup(self.tmp.cleanup)
    self.options = make.Options()
    self.options.mdir = os.path.join(self.tmp.name, "model directory")
    self.options.min_points = 1
    self.options.nthreads = 1
    with flagsaver.flagsaver(n_estimators=3, max_depth=2, nthreads=1):
      make.apply_parameter_flags(self.options)

  def tables(self, names=None, extra_activity=False, wide=False):
    names = names or [f"{i:03}" for i in range(10)]
    desc = pd.DataFrame({"Id": names, "f": np.arange(len(names), dtype=float)})
    if wide:
      desc = pd.DataFrame({"Id": names, **{
          f"f{i}": np.arange(len(names), dtype=int if i % 2 else float)
          for i in range(400)}})
    activity = pd.DataFrame({"Id": names, "response": np.arange(len(names), dtype=float)})
    if extra_activity:
      activity["other"] = 0
    dp, ap = Path(self.tmp.name)/"desc.txt", Path(self.tmp.name)/"activity.txt"
    desc.to_csv(dp, sep=" ", index=False)
    activity.to_csv(ap, sep=" ", index=False)
    return str(dp), str(ap)

  def test_wide_merge_and_saved_model(self):
    dp, ap = self.tables(wide=True)
    with warnings.catch_warnings():
      warnings.simplefilter("error", pd.errors.PerformanceWarning)
      self.assertTrue(make.build_xgboost_model(dp, ap, self.options))
    self.assertTrue((Path(self.options.mdir)/"xgboost.json").exists())
    self.assertIn("000", (Path(self.options.mdir)/"train.xy").read_text())

  def test_duplicate_and_missing_ids(self):
    dp, ap = self.tables(names=["a", "a", "b"])
    with self.assertRaises(pd.errors.MergeError):
      make.build_xgboost_model(dp, ap, self.options)
    dp, ap = self.tables()
    activity = pd.read_csv(ap, sep=" ", dtype={"Id": str})
    activity.iloc[:-1].to_csv(ap, sep=" ", index=False)
    self.assertFalse(make.build_xgboost_model(dp, ap, self.options))

  def test_extra_activity_and_too_few_rows(self):
    dp, ap = self.tables(extra_activity=True)
    self.assertFalse(make.build_xgboost_model(dp, ap, self.options))
    dp, ap = self.tables()
    self.options.min_points = 20
    self.assertFalse(make.build_xgboost_model(dp, ap, self.options))
    self.assertFalse(Path(self.options.mdir).exists())

  def test_label_mapping_matches_training(self):
    os.makedirs(self.options.mdir)
    for labels in [["A", "B", "B"], ["A", "A", "B"], [1, 2, 2]]:
      with self.subTest(labels=labels), mock.patch.object(make, "XGBClassifier") as classifier:
        self.assertTrue(make.classification(pd.DataFrame({"f": [0, 1, 2]}),
                                            np.array(labels), self.options))
        fitted = classifier.return_value.fit.call_args.args[1]
        proto = class_label_translation_pb2.ClassLabelTranslation()
        proto.ParseFromString((Path(self.options.mdir)/"class_label_translation.dat").read_bytes())
        self.assertEqual(list(fitted), [proto.to_numeric[str(v)] for v in labels])
        self.assertEqual(sum(proto.class_count.values()), 3)
        self.assertEqual(classifier.call_args.kwargs["n_estimators"], 3)

  def test_explicit_flags_override_proto(self):
    self.options.proto.eta = 0.1
    self.options.proto.tree_method = xgboost_model_pb2.EXACT
    with flagsaver.flagsaver():
      make.FLAGS(["test", "--eta=0.2", "--tree_method=hist", "--min_child_weight=7"])
      make.apply_parameter_flags(self.options)
      self.assertAlmostEqual(self.options.proto.eta, 0.2)
      self.assertEqual(make.model_parameters(self.options)["min_child_weight"], 7)
      self.assertEqual(self.options.proto.tree_method, xgboost_model_pb2.HIST)

  def test_optuna_refits_and_updates_metadata(self):
    import optuna
    real_create = optuna.create_study
    # Exercise actual CV and final XGBoost fitting with one small trial.
    def create(**kwargs):
      study = real_create(**kwargs)
      optimize = study.optimize
      study.optimize = lambda objective, **unused: optimize(objective, n_trials=1)
      study.enqueue_trial({"learning_rate": .3, "max_depth": 2,
                           "min_child_weight": 1., "subsample": 1.,
                           "colsample_bytree": 1., "gamma": 0.,
                           "reg_alpha": 1e-10, "reg_lambda": 1.})
      return study
    self.options.optuna = True
    dp, ap = self.tables()
    with mock.patch.object(optuna, "create_study", side_effect=create):
      self.assertTrue(make.build_xgboost_model(dp, ap, self.options))
    model = make.XGBRegressor()
    model.load_model(Path(self.options.mdir)/"xgboost.json")
    self.assertEqual(model.get_booster().num_boosted_rounds(), self.options.proto.n_estimators)
    proto = xgboost_model_pb2.XGBoostModel()
    proto.ParseFromString((Path(self.options.mdir)/"model_metadata.dat").read_bytes())
    self.assertAlmostEqual(proto.parameters.eta, .3)
    self.assertEqual(proto.parameters.n_estimators, model.get_booster().num_boosted_rounds())

  def test_rescore_failure_and_space_arguments(self):
    os.makedirs(self.options.mdir)
    self.options.descriptor_fname = "descriptor file.txt"
    with mock.patch.object(make.subprocess, "run", side_effect=subprocess.CalledProcessError(2, "evaluate")) as run:
      with self.assertRaises(subprocess.CalledProcessError):
        make.rescore_training_set(self.options)
      self.assertEqual(run.call_args.args[0][-1], "descriptor file.txt")
      self.assertTrue(run.call_args.kwargs["check"])

  def test_saved_classification_model(self):
    dp, ap = self.tables()
    activity = pd.read_csv(ap, sep=" ", dtype={"Id": str})
    activity["response"] = ["active"] * 3 + ["inactive"] * 7
    activity.to_csv(ap, sep=" ", index=False)
    self.options.classification = True
    self.assertTrue(make.build_xgboost_model(dp, ap, self.options))
    model = make.XGBClassifier()
    model.load_model(Path(self.options.mdir)/"xgboost.json")
    proto = class_label_translation_pb2.ClassLabelTranslation()
    proto.ParseFromString((Path(self.options.mdir)/"class_label_translation.dat").read_bytes())
    labels = {number: label for label, number in proto.to_numeric.items()}
    predictions = model.predict(pd.read_csv(dp, sep=" ").iloc[:, 1:])
    self.assertEqual(len(predictions), 10)
    self.assertTrue(all(labels[int(number)] in {"active", "inactive"} for number in predictions))
    self.assertEqual(proto.to_numeric["active"], 1)

  def test_proto_aliases_are_preserved(self):
    self.options.proto.learning_rate = .15
    with flagsaver.flagsaver():
      make.FLAGS(["test"])
      make.apply_parameter_flags(self.options)
      self.assertAlmostEqual(self.options.proto.eta, .15)

  def test_cli_failure_exit_status(self):
    result = subprocess.run([sys.executable, "-m", "xgbd.xgbd_make"],
                            capture_output=True, text=True)
    self.assertEqual(result.returncode, 1, result.stderr)


if __name__ == "__main__":
  unittest.main()

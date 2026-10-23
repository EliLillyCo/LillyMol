# Written by ChatGPT
import numpy as np
import optuna
from optuna.samplers import TPESampler
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
import xgboost as xgb


def tune_xgb_regression_optuna(
    X,
    y,
    n_trials: int = 100,
    n_splits: int = 5,
    random_state: int = 42,
    timeout: int | None = None,
):
    """
    Hyperparameter tuning for XGBoost regression using Optuna + CV + early stopping.

    Parameters
    ----------
    X, y : array-like
        Training data.
    n_trials : int
        Number of Optuna trials.
    n_splits : int
        CV folds.
    random_state : int
        RNG seed.
    timeout : int | None
        Optional time limit (seconds) for the whole study.

    Returns
    -------
    study : optuna.study.Study
        The Optuna study containing best params etc.
    """

    X = np.asarray(X)
    y = np.asarray(y)

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    def objective(trial: optuna.Trial) -> float:
        params = {
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
            "tree_method": "hist",

            # Objective / metric
            "objective": "reg:squarederror",
            "eval_metric": "rmse",

            # Repro / speed
            "random_state": random_state,
            "n_jobs": -1,
        }

        fold_rmses = []
        for fold_idx, (tr_idx, va_idx) in enumerate(kf.split(X), start=1):
            X_tr, X_va = X[tr_idx], X[va_idx]
            y_tr, y_va = y[tr_idx], y[va_idx]

            model = xgb.XGBRegressor(**params)

            # Early stopping: use a validation set within each CV fold
            model.fit(
                X_tr,
                y_tr,
                eval_set=[(X_va, y_va)],
                verbose=False,
                early_stopping_rounds=200,
            )

            preds = model.predict(X_va)
            rmse = mean_squared_error(y_va, preds, squared=False)
            fold_rmses.append(rmse)

            # Let Optuna prune unpromising trials
            trial.report(float(np.mean(fold_rmses)), step=fold_idx)
            if trial.should_prune():
                raise optuna.TrialPruned()

        return float(np.mean(fold_rmses))

    study = optuna.create_study(
        direction="minimize",
        sampler=TPESampler(seed=random_state),
        pruner=optuna.pruners.MedianPruner(n_startup_trials=10, n_warmup_steps=1),
    )

    study.optimize(objective, n_trials=n_trials, timeout=timeout, show_progress_bar=True)

    return study


# -------------------------
# Example usage
# -------------------------
# study = tune_xgb_regression_optuna(X_train, y_train, n_trials=100, n_splits=5)
# print("Best RMSE:", study.best_value)
# print("Best params:", study.best_params)

# Fit final model on full training data with best params
def fit_final_xgb_regressor(X_train, y_train, best_params: dict, random_state: int = 42):
    params = dict(best_params)
    params.update({
        "n_estimators": 10_000,
        "objective": "reg:squarederror",
        "eval_metric": "rmse",
        "tree_method": "hist",
        "random_state": random_state,
        "n_jobs": -1,
    })
    model = xgb.XGBRegressor(**params)
    # If you have a dedicated validation set, pass it here for early stopping; otherwise fit without it.
    model.fit(X_train, y_train, verbose=False)
    return model

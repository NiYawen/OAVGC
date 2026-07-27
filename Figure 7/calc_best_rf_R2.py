#!/usr/bin/env python3

import argparse
import logging
import numpy as np
import pandas as pd
import pickle

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
from joblib import Parallel, delayed


# ======================================================
# argparse
# ======================================================
def get_args():
    parser = argparse.ArgumentParser(
        description="RF forward feature selection for continuous outcome (repeated pooled R2)"
    )

    parser.add_argument("-X", required=True,
                        help="Feature matrix (rows=features, cols=samples)")
    parser.add_argument("-Y", required=True,
                        help="Metadata (rows=samples, cols=variables)")
    parser.add_argument("-y", required=True,
                        help="Continuous target column name")
    parser.add_argument("-o", required=True,
                        help="Output prefix")
    
    parser.add_argument("--fold", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=10)
    
    parser.add_argument("--n_estimators", type=int, default=500)
    parser.add_argument("--max_depth", type=int, default=10)
    parser.add_argument("--min_samples_leaf", type=int, default=3)
    
    parser.add_argument("-p", "--threads", type=int, default=-1)
    parser.add_argument("--seed", type=int, default=2026)
    
    parser.add_argument("--max_features", type=int, default=100)
    parser.add_argument("--max_corr", type=float, default=0.6)
    parser.add_argument("--step", type=float, default=1.01)
    parser.add_argument("--step_count", type=int, default=10)
    
    parser.add_argument("--min_abundance", type=float, default=0)
    parser.add_argument("--min_prevalence", type=int, default=0)
    
    parser.add_argument("--test_X", default=None,
                        help="Optional test feature matrix (rows=features, cols=samples)")
    parser.add_argument("--test_Y", default=None,
                        help="Optional test metadata (rows=samples, must contain -y  column)")
    
    parser.add_argument(
        "-s", type=int, choices=[0, 1, 2, 3], default=1,
        help="0: whitespace | 1: tab | 2: | | 3: comma"
    )
    
    return parser.parse_args()


# ======================================================
# utility
# ======================================================
def align_xy(df, y, min_abundance=0.01, min_prevalence=5):

    y = y.dropna()
    intersect = df.columns.intersection(y.index)
    
    if len(intersect) == 0:
        raise ValueError("No common samples between X and y")
    
    df = df.loc[:, intersect]
    y = y.loc[intersect]
    
    prevalence = (df > min_abundance).sum(axis=1)
    df = df.loc[prevalence >= min_prevalence]
    
    X = df.T
    
    logging.info(
        f"[ALIGN] samples={X.shape[0]}, features={X.shape[1]}"
    )
    
    return X, y


def add_feature(selected, candidate, X, cutoff=0.6):
    for f in selected:
        r, _ = spearmanr(X[f], X[candidate])
        if not np.isnan(r) and abs(r) >= cutoff:
            return False
    selected.append(candidate)
    return True


# ======================================================
# single repeat (for parallel)
# ======================================================
def _one_repeat_r2(r, X, y, n_splits, rf_params, base_seed):
    params = rf_params.copy()
    params.pop("base_seed", None)
    params["random_state"] = base_seed + r
    params["n_jobs"] = 1  # RF 内部单线程

    kf = KFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=base_seed + r
    )
    
    pred_sum = np.zeros(len(y))
    pred_n = np.zeros(len(y))
    
    for train, test in kf.split(X):
        model = RandomForestRegressor(**params)
        model.fit(X.iloc[train], y.iloc[train])
        pred = model.predict(X.iloc[test])
        pred_sum[test] += pred
        pred_n[test] += 1
    
    mean_pred = pred_sum / pred_n
    r2 = r2_score(y, mean_pred)
    
    return r + 1, r2, mean_pred


# ======================================================
# repeated pooled R2 with parallel repeats
# ======================================================
def repeated_pooled_r2_parallel(X, y, n_splits, n_repeats, rf_params):
    X = X.reset_index(drop=True)
    y = y.reset_index(drop=True)

    base_seed = rf_params.get("base_seed", 2026)
    threads = rf_params.get("n_jobs", 1)
    
    actual_threads = min(
        threads if threads > 0 else n_repeats,
        n_repeats
    )
    
    results = Parallel(n_jobs=actual_threads)(
        delayed(_one_repeat_r2)(r, X, y, n_splits, rf_params, base_seed)
        for r in range(n_repeats)
    )
    
    r2_records = [{"repeat": r, "R2": r2} for r, r2, _ in results]
    pred_all = [p for _, _, p in results]
    
    return pd.DataFrame(r2_records), np.vstack(pred_all).T

def overall_importance(X, y, args):

    model = RandomForestRegressor(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
        random_state=args.seed,
        n_jobs=args.threads
    )
    
    model.fit(X, y)
    
    return (
        pd.DataFrame({
            "feature": X.columns,
            "importance": model.feature_importances_
        })
        .sort_values("importance", ascending=False)
    )


# ======================================================
# main
# ======================================================
def main(args):

    sep_map = {0: r"\s+", 1: "\t", 2: "|", 3: ","}
    
    logging.info("[START] load data")
    df = pd.read_csv(args.X, sep=sep_map[args.s], index_col=0)
    meta = pd.read_csv(args.Y, sep=sep_map[args.s], index_col=0)
    
    y = meta[args.y].astype(float)
    
    X, y = align_xy(
        df, y,
        args.min_abundance,
        args.min_prevalence
    )
    
    logging.info(
        f"[DATA] X={X.shape[0]} samples × {X.shape[1]} features"
    )
    
    # ---------- overall importance ----------
    logging.info("[STEP] calculate overall importance")
    imp_df = overall_importance(X, y, args)
    imp_df.to_csv(f"{args.o}.overall_importance.csv", index=False)
    
    features = imp_df["feature"].tolist()
    logging.info(f"[INFO] ranked {len(features)} candidate features")
    
    # ---------- RF params ----------
    rf_params = dict(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
        n_jobs=args.threads,
        base_seed=args.seed
    )
    
    # ---------- init ----------
    modules = [features[0]]
    logging.info(f"[INIT] start with feature: {modules[0]}")
    
    r2_df, _ = repeated_pooled_r2_parallel(
        X[modules], y,
        args.fold, args.repeats,
        rf_params.copy()
    )
    
    last_r2 = r2_df["R2"].mean()
    
    logging.info(
        f"[INIT] initial pooled R2 = {last_r2:.4f} "
        f"(fold={args.fold}, repeats={args.repeats})"
    )
    
    r2_records = [
        {"feature": modules[0], "repeat": r, "R2": v, "n_features": 1}
        for r, v in zip(r2_df["repeat"], r2_df["R2"])
    ]
    
    # ---------- forward selection ----------
    fail = 0
    total = len(features) - 1
    
    for i, f in enumerate(features[1:], start=1):
    
        if len(modules) >= args.max_features or fail >= args.step_count:
            logging.info(
                f"[STOP] modules={len(modules)}, fail={fail}"
            )
            break
    
        logging.info(
            f"[TRY {i}/{total}] feature={f} | "
            f"modules={len(modules)} | fail={fail}"
        )
    
        if not add_feature(modules, f, X, args.max_corr):
            logging.info(
                f"[SKIP] feature={f} | corr > {args.max_corr}"
            )
            continue
    
        r2_new, _ = repeated_pooled_r2_parallel(
            X[modules], y,
            args.fold, args.repeats,
            rf_params.copy()
        )
    
        mean_r2 = r2_new["R2"].mean()
        threshold = last_r2 * args.step
    
        if mean_r2 >= threshold:
            logging.info(
                f"[ACCEPT] feature={f} | "
                f"R2={mean_r2:.4f} ≥ {threshold:.4f}"
            )
    
            last_r2 = mean_r2
            fail = 0
    
            for r, v in zip(r2_new["repeat"], r2_new["R2"]):
                r2_records.append(
                    {"feature": f, "repeat": r, "R2": v, "n_features": len(modules)}
                )
        else:
            modules.pop()
            fail += 1
    
            logging.info(
                f"[REJECT] feature={f} | "
                f"R2={mean_r2:.4f} < {threshold:.4f}"
            )
    
    # ---------- output ----------
    pd.DataFrame(r2_records).to_csv(
        f"{args.o}.r2_per_repeat.csv", index=False
    )
    
    final_r2, pred_final = repeated_pooled_r2_parallel(
        X[modules], y,
        args.fold, args.repeats,
        rf_params.copy()
    )
    
    final_pred = pd.DataFrame({
        "Sample": X.index,
        "mean_pred": pred_final.mean(axis=1),
        "label": y.values
    })
    final_pred.to_csv(f"{args.o}.best_pred.csv", index=False)
    
    r2_df_all = pd.DataFrame(r2_records)
    
    r2_final_df = (
        r2_df_all
        .groupby("feature", as_index=False)
        .agg(
            mean_r2=("R2", "mean"),
            n_features=("n_features", "first")
        )
    )
    
    r2_final_df.to_csv(
        f"{args.o}.r2_final.csv",
        index=False
    )
    
    logging.info(
        f"[SUMMARY] wrote R2 summary for {r2_final_df.shape[0]} steps"
    )
    
    # ---------- save final RF model ----------
    model = RandomForestRegressor(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
        random_state=args.seed,
        n_jobs=args.threads
    )
    
    model.fit(X[modules], y)
    model.myfeatures = modules
    
    with open(f"{args.o}.model", "wb") as f:
        pickle.dump(model, f)
    
    logging.info(
        f"[DONE] selected {len(modules)} features | "
        f"final pooled R2 = {final_r2['R2'].mean():.4f}"
    )
    
    # ---------- Test data ----------
    if args.test_X is not None:
        logging.info("[TEST] load test data")
    
        df_test = pd.read_csv(
            args.test_X,
            sep=sep_map[args.s],
            index_col=0
        )
    
        for f in modules:
            if f not in df_test.index:
                df_test.loc[f] = 0
    
        X_test = df_test.loc[modules].T
    
        if args.test_Y is not None:
            meta_test = pd.read_csv(
                args.test_Y,
                sep=sep_map[args.s],
                index_col=0,
                keep_default_na=False
            )
            y_test = meta_test[args.y].replace("", np.nan).astype(float)
            y_test = y_test.reindex(X_test.index)
        else:
            y_test = pd.Series(index=X_test.index, dtype=float)
    
        pred_value = model.predict(X_test)
    
        test_pred_df = pd.DataFrame({
            "Sample": X_test.index,
            "mean_pred": pred_value,
            "label": y_test.values
        })
    
        test_pred_df.to_csv(
            f"{args.o}.test_pred.csv",
            index=False
        )
    
        logging.info(
            f"[TEST] predicted {X_test.shape[0]} samples "
            f"using {len(modules)} selected features"
        )

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    args = get_args()
    main(args)

"""Phase 5 — 5-Fold Cross-Validation (การทดลอง E1–E5)
E1 random split | E2 scaffold split | E3/E4 ablation | E5 leakage demo
"""
import time

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, cross_validate, train_test_split
from sklearn.preprocessing import StandardScaler

from config import (CV_FOLDS, DATA_PROC, DESC_COLS, FEATURE_COLS, FP_COLS,
                    RANDOM_STATE, TABLES, TEST_SIZE)
from pipelines import get_estimators, make_pipe, scaffold_split

SCORING = {"r2": "r2", "rmse": "neg_root_mean_squared_error",
           "mae": "neg_mean_absolute_error"}
CV = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=RANDOM_STATE)

df = pd.read_parquet(DATA_PROC / "features_hybrid.parquet")
y = df["pIC50"].values

# ---------- splits ----------
idx_tr, idx_te = train_test_split(
    np.arange(len(df)), test_size=TEST_SIZE, random_state=RANDOM_STATE,
    stratify=df["bioactivity_class"])
sc_tr, sc_te = scaffold_split(df["canonical_smiles"].tolist(), TEST_SIZE)

np.save(DATA_PROC / "split_random.npy", np.array([idx_tr, idx_te], dtype=object),
        allow_pickle=True)
np.save(DATA_PROC / "split_scaffold.npy", np.array([sc_tr, sc_te], dtype=object),
        allow_pickle=True)

print(f"random split   : train={len(idx_tr)}  test={len(idx_te)}")
print(f"scaffold split : train={len(sc_tr)}  test={len(sc_te)}")


def run_cv(X, yv, estimators, fp_cols, desc_cols, exp, split):
    rows = []
    for name, est in estimators.items():
        t0 = time.time()
        pipe = make_pipe(est, fp_cols=fp_cols, desc_cols=desc_cols)
        cv = cross_validate(pipe, X, yv, cv=CV, scoring=SCORING,
                            return_train_score=True, n_jobs=1)
        rows.append({
            "experiment": exp, "split": split, "model": name,
            "cv_r2": cv["test_r2"].mean(), "cv_r2_sd": cv["test_r2"].std(),
            "cv_rmse": -cv["test_rmse"].mean(), "cv_rmse_sd": cv["test_rmse"].std(),
            "cv_mae": -cv["test_mae"].mean(),
            "train_r2": cv["train_r2"].mean(),
            "gap_train_cv": cv["train_r2"].mean() - cv["test_r2"].mean(),
            "fit_time_s": round(time.time() - t0, 1),
        })
        r = rows[-1]
        print(f"  {name:<17} CV R2={r['cv_r2']:.4f}±{r['cv_r2_sd']:.4f}  "
              f"RMSE={r['cv_rmse']:.4f}  MAE={r['cv_mae']:.4f}  "
              f"gap={r['gap_train_cv']:+.3f}  [{r['fit_time_s']}s]")
    return rows


results = []

print("\n=== E1: Random split, 1028 hybrid features ===")
results += run_cv(df.iloc[idx_tr], y[idx_tr], get_estimators(),
                  FP_COLS, DESC_COLS, "E1", "random")

print("\n=== E2: Scaffold split, 1028 hybrid features ===")
results += run_cv(df.iloc[sc_tr], y[sc_tr], get_estimators(),
                  FP_COLS, DESC_COLS, "E2", "scaffold")

print("\n=== E3 (ablation): Morgan FP 1024 only, Random Forest ===")
results += run_cv(df.iloc[idx_tr], y[idx_tr],
                  {"Random Forest": get_estimators()["Random Forest"]},
                  FP_COLS, [], "E3", "random")

print("\n=== E4 (ablation): Lipinski 4 descriptors only, Random Forest ===")
results += run_cv(df.iloc[idx_tr], y[idx_tr],
                  {"Random Forest": get_estimators()["Random Forest"]},
                  [], DESC_COLS, "E4", "random")

# ---------- E5: leakage demonstration ----------
print("\n=== E5: Data Leakage demonstration (O2) ===")
X_all = df[FEATURE_COLS].copy()
leaked = X_all.copy()
leaked[DESC_COLS] = StandardScaler().fit_transform(X_all[DESC_COLS])  # ผิด: fit ทั้งชุด
rows_e5 = []
for label, data, use_pipe in [("A correct (scale inside Pipeline)", df, True),
                              ("B wrong (scale before split)", leaked, False)]:
    est = get_estimators()["SVR (RBF)"]
    if use_pipe:
        pipe = make_pipe(est, FP_COLS, DESC_COLS)
        Xs = data.iloc[idx_tr]
    else:
        from sklearn.pipeline import Pipeline
        pipe = Pipeline([("model", est)])
        Xs = data.iloc[idx_tr]
    cv = cross_validate(pipe, Xs, y[idx_tr], cv=CV, scoring=SCORING, n_jobs=1)
    rows_e5.append({"experiment": "E5", "split": "random", "model": label,
                    "cv_r2": cv["test_r2"].mean(), "cv_r2_sd": cv["test_r2"].std(),
                    "cv_rmse": -cv["test_rmse"].mean(),
                    "cv_mae": -cv["test_mae"].mean()})
    print(f"  {label:<36} CV R2={cv['test_r2'].mean():.4f}  "
          f"RMSE={-cv['test_rmse'].mean():.4f}")
results += rows_e5

# ---------- E5b: feature-selection leakage (แบบเดียวกับที่พบใน xfinal.csv) ----------
print("\n=== E5b: Feature-selection leakage — เลียนแบบ xfinal.csv (k=300) ===")
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import Pipeline as SkPipeline

K = 300
rf = get_estimators()["Random Forest"]

# กรณี A (ถูก): selection อยู่ใน Pipeline -> fit เฉพาะ train fold
pipe_a = SkPipeline([("prep", make_pipe(rf, FP_COLS, DESC_COLS).named_steps["prep"]),
                     ("sel", SelectKBest(f_regression, k=K)),
                     ("model", rf)])
cv_a = cross_validate(pipe_a, df.iloc[idx_tr], y[idx_tr], cv=CV,
                      scoring=SCORING, n_jobs=1)

# กรณี B (ผิด): selection ทำบนข้อมูลทั้งชุดก่อน CV
sel_full = SelectKBest(f_regression, k=K).fit(df[FEATURE_COLS], y)
cols_leaked = np.array(FEATURE_COLS)[sel_full.get_support()]
pipe_b = SkPipeline([("model", rf)])
cv_b = cross_validate(pipe_b, df.iloc[idx_tr][cols_leaked], y[idx_tr], cv=CV,
                      scoring=SCORING, n_jobs=1)

for label, cv in [("A correct (SelectKBest inside Pipeline)", cv_a),
                  ("B wrong (SelectKBest on full dataset)", cv_b)]:
    results.append({"experiment": "E5b", "split": "random", "model": label,
                    "cv_r2": cv["test_r2"].mean(), "cv_r2_sd": cv["test_r2"].std(),
                    "cv_rmse": -cv["test_rmse"].mean(),
                    "cv_mae": -cv["test_mae"].mean()})
    print(f"  {label:<42} CV R2={cv['test_r2'].mean():.4f}  "
          f"RMSE={-cv['test_rmse'].mean():.4f}")
print(f"  -> optimistic bias = {cv_b['test_r2'].mean() - cv_a['test_r2'].mean():+.4f} R2")

out = pd.DataFrame(results)
out.to_csv(TABLES / "tbl02_cv_results.csv", index=False)
print(f"\nsaved {TABLES / 'tbl02_cv_results.csv'}")

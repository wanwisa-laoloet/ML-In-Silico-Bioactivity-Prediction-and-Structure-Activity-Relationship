"""Phase 6 — Hyperparameter Tuning และการประเมินขั้นสุดท้ายบน Test Set
Test Set ถูกแตะเพียงครั้งเดียวหลัง tuning เสร็จ
"""
import json
import time

import joblib
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, r2_score, root_mean_squared_error
from sklearn.model_selection import KFold, RandomizedSearchCV

from config import (CV_FOLDS, DATA_PROC, DESC_COLS, FIGURES, FP_COLS, MODELS,
                    N_JOBS, RANDOM_STATE, TABLES)
from pipelines import get_estimators, make_pipe

sns.set_theme(style="whitegrid", context="talk")
CV = KFold(n_splits=CV_FOLDS, shuffle=True, random_state=RANDOM_STATE)

df = pd.read_parquet(DATA_PROC / "features_hybrid.parquet")
y = df["pIC50"].values
idx_tr, idx_te = np.load(DATA_PROC / "split_random.npy", allow_pickle=True)
sc_tr, sc_te = np.load(DATA_PROC / "split_scaffold.npy", allow_pickle=True)

GRIDS = {
    "Random Forest": {
        "model__n_estimators": [300, 500, 800],
        "model__max_depth": [None, 15, 25],
        "model__min_samples_split": [2, 5, 10],
        "model__min_samples_leaf": [1, 2, 4],
        "model__max_features": ["sqrt", 0.3],
    },
    "XGBoost": {
        "model__n_estimators": [300, 600],
        "model__learning_rate": [0.03, 0.05, 0.1],
        "model__max_depth": [4, 6, 8],
        "model__subsample": [0.8, 1.0],
        "model__colsample_bytree": [0.6, 0.8],
    },
    "SVR (RBF)": {
        "model__C": [1, 10, 100],
        "model__gamma": ["scale", 0.001, 0.01],
        "model__epsilon": [0.05, 0.1, 0.2],
    },
}
N_ITER = {"Random Forest": 30, "XGBoost": 30, "SVR (RBF)": 27}


def metrics(yt, yp):
    return {"r2": r2_score(yt, yp), "rmse": root_mean_squared_error(yt, yp),
            "mae": mean_absolute_error(yt, yp),
            "spearman": spearmanr(yt, yp).statistic}


def evaluate(split_name, tr, te):
    print(f"\n{'=' * 62}\nTUNING on {split_name} split "
          f"(train={len(tr)}, test={len(te)})\n{'=' * 62}")
    rows, best_models, preds = [], {}, {}

    # baseline: Dummy (ไม่ต้อง tune)
    dummy = make_pipe(get_estimators()["Dummy Regressor"])
    dummy.fit(df.iloc[tr], y[tr])
    yp = dummy.predict(df.iloc[te])
    rows.append({"split": split_name, "model": "Dummy Regressor",
                 "best_params": "-", "fit_time_s": 0.0, **metrics(y[te], yp)})
    preds["Dummy Regressor"] = yp
    print(f"  Dummy Regressor   Test R2={rows[-1]['r2']:.4f}")

    for name, grid in GRIDS.items():
        t0 = time.time()
        search = RandomizedSearchCV(
            make_pipe(get_estimators()[name]), grid, n_iter=N_ITER[name],
            cv=CV, scoring="r2", random_state=RANDOM_STATE, n_jobs=1, refit=True)
        search.fit(df.iloc[tr], y[tr])
        yp = search.best_estimator_.predict(df.iloc[te])
        m = metrics(y[te], yp)
        el = round(time.time() - t0, 1)
        rows.append({"split": split_name, "model": name,
                     "best_cv_r2": search.best_score_,
                     "best_params": json.dumps(
                         {k.replace('model__', ''): v
                          for k, v in search.best_params_.items()}),
                     "fit_time_s": el, **m})
        best_models[name] = search.best_estimator_
        preds[name] = yp
        print(f"  {name:<17} Test R2={m['r2']:.4f}  RMSE={m['rmse']:.4f}  "
              f"MAE={m['mae']:.4f}  rho={m['spearman']:.4f}  [{el}s]")
        print(f"      best: {rows[-1]['best_params']}")
    return rows, best_models, preds


rows_r, models_r, preds_r = evaluate("random", idx_tr, idx_te)
rows_s, models_s, preds_s = evaluate("scaffold", sc_tr, sc_te)

res = pd.DataFrame(rows_r + rows_s)
res.to_csv(TABLES / "tbl03_test_results.csv", index=False)

# ---------- เลือกและบันทึกโมเดลที่ดีที่สุด (จาก random split ตามข้อกำหนดใน PDF) ----------
tuned = res[(res["split"] == "random") & (res["model"] != "Dummy Regressor")]
best_name = tuned.loc[tuned["r2"].idxmax(), "model"]
best_model = models_r[best_name]
stamp = time.strftime("%Y%m%d")
slug = best_name.split()[0].lower()
joblib.dump(best_model, MODELS / f"best_model_{slug}_{stamp}.joblib")
print(f"\n>>> BEST MODEL: {best_name}  (saved as best_model_{slug}_{stamp}.joblib)")

pd.DataFrame({
    "canonical_smiles": df.iloc[idx_te]["canonical_smiles"].values,
    "inchikey": df.iloc[idx_te]["inchikey"].values,
    "y_true": y[idx_te], "y_pred": preds_r[best_name],
    "residual": y[idx_te] - preds_r[best_name],
    "bioactivity_class": df.iloc[idx_te]["bioactivity_class"].values,
}).to_csv(DATA_PROC / "predictions_test.csv", index=False)

# ---------- fig03: model comparison ----------
cv_tbl = pd.read_csv(TABLES / "tbl02_cv_results.csv")
e1 = cv_tbl[cv_tbl["experiment"] == "E1"].set_index("model")
order = ["Dummy Regressor", "Random Forest", "XGBoost", "SVR (RBF)"]
rr = res[res["split"] == "random"].set_index("model")

fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))
for ax, (metric, cvk, label) in zip(axes, [("r2", "cv_r2", "R²"),
                                           ("rmse", "cv_rmse", "RMSE"),
                                           ("mae", "cv_mae", "MAE")]):
    x = np.arange(len(order))
    ax.bar(x - .2, [e1.loc[m, cvk] for m in order], .4, label="5-Fold CV",
           color="#8ecae6", edgecolor="#023047")
    ax.bar(x + .2, [rr.loc[m, metric] for m in order], .4, label="Test (25%)",
           color="#fb8500", edgecolor="#023047")
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace(" ", "\n") for m in order], fontsize=11)
    ax.set_title(label)
    if metric == "r2":
        ax.axhline(0, color="black", lw=1)
        ax.legend(fontsize=11)
fig.suptitle("Model comparison — CV vs held-out Test set (random split)", y=1.03)
fig.tight_layout()
fig.savefig(FIGURES / "fig03_model_comparison.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------- fig04: residual analysis ----------
yt, yp = y[idx_te], preds_r[best_name]
resid = yt - yp
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
axes[0].scatter(yt, yp, s=26, alpha=.55, color="#264653", edgecolors="none")
lo, hi = min(yt.min(), yp.min()) - .3, max(yt.max(), yp.max()) + .3
axes[0].plot([lo, hi], [lo, hi], "--", color="#e76f51", lw=2, label="y = x")
axes[0].set(xlabel="Actual pIC50", ylabel="Predicted pIC50",
            title=f"{best_name} — Predicted vs Actual\n"
                  f"R² = {rr.loc[best_name, 'r2']:.3f}, "
                  f"RMSE = {rr.loc[best_name, 'rmse']:.3f}")
axes[0].legend(fontsize=11)
axes[1].scatter(yp, resid, s=26, alpha=.55, color="#264653", edgecolors="none")
axes[1].axhline(0, color="#e76f51", ls="--", lw=2)
axes[1].axhline(1, color="grey", ls=":", lw=1.5)
axes[1].axhline(-1, color="grey", ls=":", lw=1.5)
axes[1].set(xlabel="Predicted pIC50", ylabel="Residual (actual − predicted)",
            title=f"Residual plot — {(np.abs(resid) <= 1).mean() * 100:.1f}% "
                  f"within ±1 log unit")
fig.tight_layout()
fig.savefig(FIGURES / "fig04_residual_plot.png", dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 5.5))
sns.histplot(resid, bins=40, kde=True, color="#264653", ax=ax)
ax.axvline(0, color="#e76f51", ls="--", lw=2)
ax.set(xlabel="Residual", ylabel="Count",
       title=f"Residual distribution (mean={resid.mean():.3f}, sd={resid.std():.3f})")
fig.tight_layout()
fig.savefig(FIGURES / "fig04b_residual_hist.png", dpi=150)
plt.close(fig)

print("\n=== Test set results ===")
print(res[["split", "model", "r2", "rmse", "mae", "spearman"]].round(4).to_string(index=False))

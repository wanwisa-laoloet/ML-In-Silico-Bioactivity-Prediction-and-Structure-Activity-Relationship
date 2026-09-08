"""Phase 8 — External Validation, BLINDED prediction และ Applicability Domain"""
import glob
import matplotlib
matplotlib.use("Agg")

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, r2_score, root_mean_squared_error

from config import (DATA_INTERIM, DATA_PROC, FEATURE_COLS, FIGURES, FP_BITS,
                    FP_RADIUS, MODELS, TABLES)
from features import build

RDLogger.DisableLog("rdApp.*")
sns.set_theme(style="whitegrid", context="talk")

model_path = sorted(glob.glob(str(MODELS / "best_model_*.joblib")))[-1]
model = joblib.load(model_path)
print(f"loaded model: {model_path.split(chr(92))[-1]}")

train = pd.read_parquet(DATA_PROC / "features_hybrid.parquet")
idx_tr, _ = np.load(DATA_PROC / "split_random.npy", allow_pickle=True)

# ---------- T8.1 external validation on DDH ----------
ddh = pd.read_parquet(DATA_PROC / "external_ddh_features.parquet")
yp = model.predict(ddh)
yt = ddh["pIC50"].values
ext = {"n": len(ddh), "r2": r2_score(yt, yp),
       "rmse": root_mean_squared_error(yt, yp), "mae": mean_absolute_error(yt, yp),
       "spearman": spearmanr(yt, yp).statistic,
       "bias_mean_error": float(np.mean(yp - yt))}
print("\n=== T8.1 DDH external validation ===")
for k, v in ext.items():
    print(f"  {k:<16} {v:.4f}" if isinstance(v, float) else f"  {k:<16} {v}")

pd.DataFrame({"compound": ddh["molecule_chembl_id"], "smiles": ddh["canonical_smiles"],
              "y_true": yt, "y_pred": yp, "residual": yt - yp}).to_csv(
    DATA_PROC / "predictions_external_ddh.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(15, 6))
axes[0].scatter(yt, yp, s=45, alpha=.7, color="#e76f51", edgecolors="none")
lo, hi = 2.8, 8.0
axes[0].plot([lo, hi], [lo, hi], "--", color="#264653", lw=2, label="y = x")
axes[0].set(xlabel="Actual pIC50 (DDH)", ylabel="Predicted pIC50", xlim=(lo, hi),
            ylim=(lo, hi), title=f"External validation — DDH (n={len(ddh)})\n"
                                 f"R² = {ext['r2']:.3f}, ρ = {ext['spearman']:.3f}")
axes[0].legend(fontsize=11)
sns.kdeplot(train["pIC50"], ax=axes[1], label=f"ChEMBL train (n={len(train)})",
            fill=True, color="#264653")
sns.kdeplot(yt, ax=axes[1], label=f"DDH external (n={len(ddh)})", fill=True,
            color="#e76f51")
axes[1].set(xlabel="pIC50", title="Label distribution shift")
axes[1].legend(fontsize=11)
fig.tight_layout()
fig.savefig(FIGURES / "fig08_external_validation.png", dpi=150)
plt.close(fig)

# ---------- T8.3 BLINDED predictions ----------
bl = pd.read_csv(DATA_INTERIM / "ddh_blinded.csv")
bl_feat, _ = build(bl, smiles_col="SMILES", collect_bitinfo=False)
bl_feat["pred_pIC50"] = model.predict(bl_feat)
bl_feat["pred_IC50_uM"] = (10 ** (9 - bl_feat["pred_pIC50"])) / 1000
out_bl = bl_feat[["compound_no", "SMILES", "pred_pIC50", "pred_IC50_uM"]].round(3)
out_bl.to_csv(TABLES / "tbl07_blinded_predictions.csv", index=False)
print(f"\n=== T8.3 BLINDED predictions (n={len(out_bl)}) ===")
print(out_bl.to_string(index=False))

# ---------- T8.4 Applicability Domain ----------
gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


def fps_of(smiles_series):
    out = []
    for s in smiles_series:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out.append(gen.GetFingerprint(m))
    return out


train_fps = fps_of(train.iloc[idx_tr]["canonical_smiles"])


def max_sim(smiles_series):
    vals = []
    for s in smiles_series:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        vals.append(max(DataStructs.BulkTanimotoSimilarity(
            gen.GetFingerprint(m), train_fps)))
    return np.array(vals)


mt = pd.read_parquet(DATA_PROC / "multitarget_features.parquet")
test_idx = np.load(DATA_PROC / "split_random.npy", allow_pickle=True)[1]
sims = {
    "ChEMBL test (in-domain)": max_sim(train.iloc[test_idx]["canonical_smiles"]),
    "DDH external": max_sim(ddh["canonical_smiles"]),
    "MultiTarget (out-of-domain)": max_sim(mt["SMILES"]),
}

rf_like = hasattr(model.named_steps["model"], "estimators_")
ad_rows = []
for name, s in sims.items():
    ad_rows.append({"dataset": name, "n": len(s),
                    "max_tanimoto_mean": round(float(s.mean()), 4),
                    "max_tanimoto_median": round(float(np.median(s)), 4),
                    "pct_below_0.4": round(100 * float((s < 0.4).mean()), 1)})
ad = pd.DataFrame(ad_rows)
ad.to_csv(TABLES / "tbl06_applicability_domain.csv", index=False)
print("\n=== T8.4 Applicability Domain (max Tanimoto to training set) ===")
print(ad.to_string(index=False))

fig, ax = plt.subplots(figsize=(10, 6))
for (name, s), c in zip(sims.items(), ["#264653", "#e76f51", "#adb5bd"]):
    sns.kdeplot(s, ax=ax, label=f"{name} (n={len(s)})", fill=True, color=c, alpha=.45)
ax.axvline(0.4, color="black", ls="--", lw=1.5)
ax.text(0.41, ax.get_ylim()[1] * .9, "AD threshold 0.4", fontsize=11)
ax.set(xlabel="Max Tanimoto similarity to training set", ylabel="Density",
       title="Applicability Domain — similarity to ChEMBL Mpro training set")
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(FIGURES / "fig09_applicability_domain.png", dpi=150)
plt.close(fig)

# ---------- prediction spread (RF only) ----------
if rf_like:
    prep = model.named_steps["prep"]
    rf = model.named_steps["model"]

    def spread(frame):
        Xp = np.asarray(prep.transform(frame))
        allp = np.stack([t.predict(Xp) for t in rf.estimators_])
        return allp.std(axis=0)

    sp = {"ChEMBL test (in-domain)": spread(train.iloc[test_idx]),
          "DDH external": spread(ddh),
          "MultiTarget (out-of-domain)": spread(mt)}
    sp_tbl = pd.DataFrame([{"dataset": k, "rf_pred_sd_mean": round(float(v.mean()), 4)}
                           for k, v in sp.items()])
    sp_tbl.to_csv(TABLES / "tbl06b_rf_prediction_spread.csv", index=False)
    print("\n=== RF prediction spread (uncertainty proxy) ===")
    print(sp_tbl.to_string(index=False))

pd.DataFrame([ext]).to_csv(TABLES / "tbl08_external_metrics.csv", index=False)
print("\nfigures:", sorted(p.name for p in FIGURES.glob("fig0[89]*")))

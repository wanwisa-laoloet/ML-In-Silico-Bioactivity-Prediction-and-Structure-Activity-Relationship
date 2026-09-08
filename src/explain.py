"""Phase 7 — SHAP และการวิเคราะห์ SAR (วัตถุประสงค์ O3)
Mean|SHAP| -> Top-20 features -> Bit Mapping -> ภาพโครงสร้างย่อย 2 มิติ
"""
import json
import pickle

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw
from sklearn.model_selection import RandomizedSearchCV, KFold

from config import (CV_FOLDS, DATA_PROC, DESC_COLS, FEATURE_COLS, FIGURES,
                    FP_COLS, RANDOM_STATE, TABLES)
from pipelines import get_estimators, make_pipe

RDLogger.DisableLog("rdApp.*")

df = pd.read_parquet(DATA_PROC / "features_hybrid.parquet")
y = df["pIC50"].values
idx_tr, idx_te = np.load(DATA_PROC / "split_random.npy", allow_pickle=True)
res = pd.read_csv(TABLES / "tbl03_test_results.csv")

# ---------- เลือกโมเดล tree ที่ดีที่สุด (TreeExplainer ต้องใช้ tree-based) ----------
tree_res = res[(res["split"] == "random")
               & (res["model"].isin(["Random Forest", "XGBoost"]))]
best_name = tree_res.loc[tree_res["r2"].idxmax(), "model"]
best_params = json.loads(tree_res.loc[tree_res["r2"].idxmax(), "best_params"])
print(f"SHAP model: {best_name}  params={best_params}")

est = get_estimators()[best_name]
est.set_params(**best_params)
pipe = make_pipe(est)
pipe.fit(df.iloc[idx_tr], y[idx_tr])

# ---------- คำนวณ SHAP บน Test Set ----------
X_test = pipe.named_steps["prep"].transform(df.iloc[idx_te])
X_test = pd.DataFrame(np.asarray(X_test), columns=FEATURE_COLS)

MAX_N = 400
rng = np.random.default_rng(RANDOM_STATE)
sub = (rng.choice(len(X_test), MAX_N, replace=False)
       if len(X_test) > MAX_N else np.arange(len(X_test)))
X_shap = X_test.iloc[sub]
print(f"computing SHAP on {len(X_shap)} test molecules ...")

explainer = shap.TreeExplainer(pipe.named_steps["model"])
sv = explainer.shap_values(X_shap, check_additivity=False)
mean_abs = np.abs(sv).mean(axis=0)  # Importance_j = (1/N) * sum|phi_j|

imp = (pd.DataFrame({"feature": FEATURE_COLS, "mean_abs_shap": mean_abs})
       .sort_values("mean_abs_shap", ascending=False).reset_index(drop=True))
imp["type"] = np.where(imp["feature"].isin(DESC_COLS), "Lipinski descriptor",
                       "Morgan fingerprint bit")
imp.head(50).to_csv(TABLES / "tbl04_shap_importance_top50.csv", index=False)

# ---------- fig05: bar chart Top-20 ----------
top20 = imp.head(20).iloc[::-1]
colors = ["#e76f51" if t == "Lipinski descriptor" else "#264653"
          for t in top20["type"]]
fig, ax = plt.subplots(figsize=(10, 8))
ax.barh(top20["feature"], top20["mean_abs_shap"], color=colors)
ax.set(xlabel="Mean |SHAP value|", ylabel="",
       title=f"Top-20 feature importance — {best_name}\n"
             f"Mean|SHAP| on {len(X_shap)} held-out test molecules")
handles = [plt.Rectangle((0, 0), 1, 1, color="#e76f51"),
           plt.Rectangle((0, 0), 1, 1, color="#264653")]
ax.legend(handles, ["Lipinski descriptor", "Morgan fingerprint bit"], fontsize=11)
fig.tight_layout()
fig.savefig(FIGURES / "fig05_shap_bar_top20.png", dpi=150)
plt.close(fig)

# ---------- fig06: beeswarm ----------
plt.figure(figsize=(10, 8))
shap.summary_plot(sv, X_shap, max_display=20, show=False)
plt.title(f"SHAP beeswarm — {best_name}", fontsize=14)
plt.tight_layout()
plt.savefig(FIGURES / "fig06_shap_beeswarm.png", dpi=150, bbox_inches="tight")
plt.close()

# ---------- fig06b: dependence plot LogP ----------
plt.figure(figsize=(9, 6))
shap.dependence_plot("LogP", sv, X_shap, interaction_index=None, show=False)
plt.title("SHAP dependence — LogP", fontsize=14)
plt.tight_layout()
plt.savefig(FIGURES / "fig06b_shap_dependence_logp.png", dpi=150,
            bbox_inches="tight")
plt.close()

# ---------- Bit Mapping -> ภาพโครงสร้างย่อย ----------
with open(DATA_PROC / "bitinfo.pkl", "rb") as f:
    bitinfo = pickle.load(f)

test_df = df.iloc[idx_te].reset_index(drop=True).iloc[sub].reset_index(drop=True)
sv_df = pd.DataFrame(sv, columns=FEATURE_COLS)

fp_imp = imp[imp["type"] == "Morgan fingerprint bit"].copy()
# ทิศทางเฉลี่ยของ SHAP เมื่อบิตถูกเซ็ต (=1)
direction = {}
for feat in fp_imp["feature"]:
    on = X_shap[feat].values == 1
    direction[feat] = sv_df.loc[on, feat].mean() if on.sum() else 0.0
fp_imp["mean_shap_when_on"] = fp_imp["feature"].map(direction)

up = fp_imp[fp_imp["mean_shap_when_on"] > 0].head(3)
down = fp_imp[fp_imp["mean_shap_when_on"] < 0].head(3)
selected = pd.concat([up, down])

sar_rows = []
for _, r in selected.iterrows():
    feat = r["feature"]
    bit_id = int(feat.split("_")[1])
    on_mask = X_shap[feat].values == 1
    if on_mask.sum() == 0:
        continue
    # โมเลกุลที่มีบิตนี้และ |SHAP| สูงสุด
    cand = np.where(on_mask)[0]
    pick = cand[np.argmax(np.abs(sv_df[feat].values[cand]))]
    smi = test_df.loc[pick, "canonical_smiles"]
    mol = Chem.MolFromSmiles(smi)
    out = FIGURES / f"fig07_bit_{bit_id}_substructure.png"
    try:
        img = Draw.DrawMorganBit(mol, bit_id, bitinfo[smi], useSVG=False)
        # RDKit รุ่นใหม่คืน PNG bytes, รุ่นเก่าคืน PIL.Image, โหมด SVG คืน str
        if isinstance(img, (bytes, bytearray)):
            out.write_bytes(img)
        elif isinstance(img, str):
            out.with_suffix(".svg").write_text(img, encoding="utf-8")
        else:
            img.save(out)
        drawn = True
    except Exception as exc:  # noqa: BLE001
        print(f"  ! cannot draw bit {bit_id}: {exc}")
        drawn = False

    full = df[feat].values == 1
    act = df["bioactivity_class"].values == "active"
    ina = df["bioactivity_class"].values == "inactive"
    sar_rows.append({
        "bit": feat, "bit_id": bit_id,
        "mean_abs_shap": round(r["mean_abs_shap"], 4),
        "direction": "increases pIC50" if r["mean_shap_when_on"] > 0
                     else "decreases pIC50",
        "mean_shap_when_on": round(r["mean_shap_when_on"], 4),
        "pct_in_active": round(100 * full[act].mean(), 1),
        "pct_in_inactive": round(100 * full[ina].mean(), 1),
        "n_molecules_with_bit": int(full.sum()),
        "example_smiles": smi,
        "figure": out.name if drawn else "",
    })

sar = pd.DataFrame(sar_rows)
sar.to_csv(TABLES / "tbl05_sar_bits.csv", index=False)

print("\n=== Top-10 features by Mean|SHAP| ===")
print(imp.head(10).to_string(index=False))
print("\n=== SAR bit table ===")
print(sar[["bit", "mean_abs_shap", "direction", "pct_in_active",
           "pct_in_inactive", "n_molecules_with_bit"]].to_string(index=False))
print("\nfigures:", sorted(p.name for p in FIGURES.glob("fig0[567]*")))

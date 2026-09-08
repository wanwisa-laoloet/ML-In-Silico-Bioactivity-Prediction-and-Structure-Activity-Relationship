"""Phase 4 — Exploratory Data Analysis
สร้างภาพ fig01–fig02c และตารางสรุปสถิติ
"""
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import mannwhitneyu, skew

from config import DATA_PROC, DESC_COLS, FIGURES, FP_COLS, TABLES

RDLogger.DisableLog("rdApp.*")
sns.set_theme(style="whitegrid", context="talk")
PALETTE = {"active": "#2a9d8f", "intermediate": "#e9c46a", "inactive": "#e76f51"}

df = pd.read_parquet(DATA_PROC / "features_hybrid.parquet")
ddh = pd.read_parquet(DATA_PROC / "external_ddh_features.parquet")
mt = pd.read_parquet(DATA_PROC / "multitarget_features.parquet")

# ---------- fig01: pIC50 distribution ----------
fig, ax = plt.subplots(figsize=(10, 6))
sns.histplot(df["pIC50"], bins=45, kde=True, color="#264653", ax=ax)
m, med, sd, sk = (df["pIC50"].mean(), df["pIC50"].median(),
                  df["pIC50"].std(), skew(df["pIC50"]))
ax.axvline(m, color="#e76f51", ls="--", lw=2, label=f"mean = {m:.2f}")
ax.axvline(med, color="#2a9d8f", ls=":", lw=2, label=f"median = {med:.2f}")
ax.set(xlabel="pIC50", ylabel="Count",
       title=f"Distribution of pIC50 — SARS-CoV-2 Mpro (n={len(df)})")
ax.legend(title=f"SD = {sd:.2f}   skew = {sk:.2f}", fontsize=11)
fig.tight_layout()
fig.savefig(FIGURES / "fig01_pic50_distribution.png", dpi=150)
plt.close(fig)

# ---------- fig02: chemical space MW vs LogP ----------
fig, ax = plt.subplots(figsize=(10, 7))
for cls in ["inactive", "intermediate", "active"]:
    s = df[df["bioactivity_class"] == cls]
    ax.scatter(s["MW"], s["LogP"], s=22, alpha=.55, label=f"{cls} (n={len(s)})",
               color=PALETTE[cls], edgecolors="none")
ax.axvline(500, color="grey", ls="--", lw=1.5)
ax.axhline(5, color="grey", ls="--", lw=1.5)
ax.text(505, 8.3, "MW = 500", color="grey", fontsize=11)
ax.text(30, 5.15, "LogP = 5", color="grey", fontsize=11)
ax.set(xlabel="Molecular Weight (Da)", ylabel="LogP (Crippen)",
       title="Chemical Space — Lipinski boundaries")
ax.legend(fontsize=11)
fig.tight_layout()
fig.savefig(FIGURES / "fig02_chemical_space_mw_logp.png", dpi=150)
plt.close(fig)

# ---------- fig02b: Lipinski descriptors by class + Mann-Whitney ----------
fig, axes = plt.subplots(1, 4, figsize=(18, 5))
stats_rows = []
order = ["inactive", "intermediate", "active"]
for ax, col in zip(axes, DESC_COLS):
    sns.boxplot(data=df, x="bioactivity_class", y=col, order=order,
                hue="bioactivity_class", palette=PALETTE, legend=False, ax=ax)
    ax.set(xlabel="", title=col)
    ax.tick_params(axis="x", rotation=20)
    a = df.loc[df["bioactivity_class"] == "active", col]
    i = df.loc[df["bioactivity_class"] == "inactive", col]
    u, p = mannwhitneyu(a, i)
    stats_rows.append({"descriptor": col, "active_median": a.median(),
                       "inactive_median": i.median(), "U": u, "p_value": p})
fig.suptitle("Lipinski descriptors by bioactivity class", y=1.02)
fig.tight_layout()
fig.savefig(FIGURES / "fig02b_lipinski_by_class.png", dpi=150, bbox_inches="tight")
plt.close(fig)
stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(TABLES / "tbl01b_lipinski_mannwhitney.csv", index=False)

# ---------- fig02c: chemical space of 3 datasets ----------
fig, ax = plt.subplots(figsize=(10, 7))
for d, name, c in [(mt, "MultiTarget (DRD2/HDAC1/BACE1)", "#adb5bd"),
                   (df, "ChEMBL Mpro (primary)", "#264653"),
                   (ddh, "DDH COVID (external)", "#e76f51")]:
    ax.scatter(d["MW"], d["LogP"], s=20, alpha=.5, label=f"{name}  n={len(d)}",
               color=c, edgecolors="none")
ax.axvline(500, color="grey", ls="--", lw=1.2)
ax.axhline(5, color="grey", ls="--", lw=1.2)
ax.set(xlabel="Molecular Weight (Da)", ylabel="LogP (Crippen)",
       title="Chemical space comparison across three datasets")
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(FIGURES / "fig02c_chemspace_3sets.png", dpi=150)
plt.close(fig)

# ---------- Murcko scaffolds ----------
scaf = df["canonical_smiles"].apply(
    lambda s: MurckoScaffold.MurckoScaffoldSmiles(smiles=s, includeChirality=False))
top = scaf.value_counts().head(10).reset_index()
top.columns = ["murcko_scaffold", "count"]
top["pct"] = (100 * top["count"] / len(df)).round(2)
top.to_csv(TABLES / "tbl01c_top_scaffolds.csv", index=False)

# ---------- summary table ----------
bits = df[FP_COLS].values
summary = pd.DataFrame([{
    "dataset": "ChEMBL Mpro (primary)", "n_molecules": len(df),
    "pIC50_mean": round(df["pIC50"].mean(), 3),
    "pIC50_sd": round(df["pIC50"].std(), 3),
    "pIC50_min": round(df["pIC50"].min(), 3),
    "pIC50_max": round(df["pIC50"].max(), 3),
    "n_active": int((df["bioactivity_class"] == "active").sum()),
    "n_intermediate": int((df["bioactivity_class"] == "intermediate").sum()),
    "n_inactive": int((df["bioactivity_class"] == "inactive").sum()),
    "unique_scaffolds": scaf.nunique(),
    "mean_bits_set": round(bits.sum(1).mean(), 1),
    "zero_variance_bits": int((df[FP_COLS].nunique() == 1).sum()),
}, {
    "dataset": "DDH COVID (external)", "n_molecules": len(ddh),
    "pIC50_mean": round(ddh["pIC50"].mean(), 3),
    "pIC50_sd": round(ddh["pIC50"].std(), 3),
    "pIC50_min": round(ddh["pIC50"].min(), 3),
    "pIC50_max": round(ddh["pIC50"].max(), 3),
    "n_active": int((ddh["bioactivity_class"] == "active").sum()),
    "n_intermediate": int((ddh["bioactivity_class"] == "intermediate").sum()),
    "n_inactive": int((ddh["bioactivity_class"] == "inactive").sum()),
    "unique_scaffolds": None, "mean_bits_set": round(ddh[FP_COLS].values.sum(1).mean(), 1),
    "zero_variance_bits": None,
}])
summary.to_csv(TABLES / "tbl01_dataset_summary.csv", index=False)

print("=== EDA summary ===")
print(summary.to_string(index=False))
print(f"\nunique Murcko scaffolds: {scaf.nunique()} / {len(df)} molecules "
      f"({100 * scaf.nunique() / len(df):.1f}%)")
print("\nTop-5 scaffolds:")
print(top.head(5).to_string(index=False))
print("\nMann-Whitney U (active vs inactive):")
print(stats_df.to_string(index=False))
print("\nfigures saved:", sorted(p.name for p in FIGURES.glob("fig0*")))

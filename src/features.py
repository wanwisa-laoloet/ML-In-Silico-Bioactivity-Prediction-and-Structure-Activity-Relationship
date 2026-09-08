"""Phase 3 — Feature Extraction
Morgan Fingerprints (r=2, 1024 bits) + Lipinski Descriptors 4 ตัว = 1,028 features
เก็บ bitInfo ไว้สำหรับ DrawMorganBit ใน Phase 7
"""
import pickle

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdFingerprintGenerator

from config import (DATA_PROC, DATASET_ORIG, DESC_COLS, FP_BITS, FP_COLS,
                    FP_RADIUS)

RDLogger.DisableLog("rdApp.*")
_gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


def featurize(smiles, keep_bitinfo=True):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitInfoMap()
    fp = _gen.GetFingerprintAsNumPy(mol, additionalOutput=ao)
    desc = [Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
            Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)]
    return fp.astype(np.uint8), desc, (ao.GetBitInfoMap() if keep_bitinfo else None)


def build(df, smiles_col="canonical_smiles", collect_bitinfo=True):
    fps, descs, bitinfo, keep = [], [], {}, []
    for i, smi in enumerate(df[smiles_col]):
        fp, desc, bi = featurize(smi, keep_bitinfo=collect_bitinfo)
        if fp is None:
            continue
        fps.append(fp)
        descs.append(desc)
        keep.append(i)
        if collect_bitinfo:
            bitinfo[smi] = {int(k): [tuple(t) for t in v] for k, v in bi.items()}

    X = pd.DataFrame(np.vstack(fps), columns=FP_COLS)
    X[DESC_COLS] = np.array(descs)
    out = df.iloc[keep].reset_index(drop=True)
    for c in X.columns:
        out[c] = X[c].values
    return out, bitinfo


if __name__ == "__main__":
    print("=== ChEMBL master ===")
    master = pd.read_csv(DATA_PROC / "clean_bioactivity_master.csv")
    feat, bitinfo = build(master)
    feat.to_parquet(DATA_PROC / "features_hybrid.parquet", index=False)
    with open(DATA_PROC / "bitinfo.pkl", "wb") as f:
        pickle.dump(bitinfo, f)
    print(f"features_hybrid.parquet shape={feat.shape}  "
          f"(features={len(FP_COLS) + len(DESC_COLS)})")

    # ตรวจคุณภาพ
    fp_vals = np.unique(feat[FP_COLS].values)
    zero_var = int((feat[FP_COLS].nunique() == 1).sum())
    print(f"  fingerprint unique values: {fp_vals}")
    print(f"  zero-variance bits: {zero_var} / {FP_BITS} "
          f"(informative: {FP_BITS - zero_var})")
    print(f"  bits set per molecule: mean={feat[FP_COLS].sum(axis=1).mean():.1f}")
    print(f"  NaN in features: {int(feat[FP_COLS + DESC_COLS].isna().sum().sum())}")
    print(f"  MW == 0 rows: {int((feat['MW'] == 0).sum())}")
    print(feat[DESC_COLS].describe().round(2).to_string())

    print("\n=== DDH external ===")
    ddh = pd.read_csv(DATA_PROC / "external_ddh_clean.csv")
    ddh_feat, _ = build(ddh, collect_bitinfo=False)
    ddh_feat.to_parquet(DATA_PROC / "external_ddh_features.parquet", index=False)
    print(f"external_ddh_features.parquet shape={ddh_feat.shape}")

    print("\n=== MultiTarget (Applicability Domain reference) ===")
    mt = pd.read_csv(DATASET_ORIG / "MultiTarget Bioactivity ChEMBL" /
                     "molecule_classification_dataset.csv")
    mt = mt.dropna(subset=["SMILES"]).drop_duplicates(subset=["SMILES"])
    mt_feat, _ = build(mt, smiles_col="SMILES", collect_bitinfo=False)
    mt_feat.to_parquet(DATA_PROC / "multitarget_features.parquet", index=False)
    print(f"multitarget_features.parquet shape={mt_feat.shape}")

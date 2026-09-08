"""Phase 2 — Data Cleaning & Standardization
ทำตามลำดับ 13 ขั้นตอนใน Project_Plan.md ส่วนที่ 6 (Phase 2) อย่างเคร่งครัด
ผลผลิต: data/processed/clean_bioactivity_master.csv
"""
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from config import (DATA_INTERIM, DATA_PROC, DEDUP_SD_MAX, HEAVY_ATOM_MAX,
                    HEAVY_ATOM_MIN, PIC50_MAX, PIC50_MIN, bioactivity_class)

RDLogger.DisableLog("rdApp.*")

_lfc = rdMolStandardize.LargestFragmentChooser()
_norm = rdMolStandardize.Normalizer()
_unch = rdMolStandardize.Uncharger()


def standardize(smiles):
    """salt strip -> normalize -> uncharge -> canonical SMILES + InChIKey"""
    if not isinstance(smiles, str) or not smiles:
        return None, None, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    try:
        mol = _unch.uncharge(_norm.normalize(_lfc.choose(mol)))
    except Exception:  # noqa: BLE001
        return None, None, None
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        return None, None, None
    return (Chem.MolToSmiles(mol, canonical=True),
            Chem.MolToInchiKey(mol),
            mol.GetNumHeavyAtoms())


class Log:
    def __init__(self):
        self.rows = []

    def step(self, name, df):
        self.rows.append({"step": name, "n_rows": len(df)})
        prev = self.rows[-2]["n_rows"] if len(self.rows) > 1 else len(df)
        print(f"  {name:<46} n={len(df):>6}  (-{prev - len(df)})")
        return df

    def frame(self):
        return pd.DataFrame(self.rows)


def clean(df, source, value_col="standard_value", pic50_col=None):
    log = Log()
    log.step("0. input", df)

    if pic50_col is None:
        df = df.dropna(subset=["canonical_smiles", value_col])
        log.step("1. drop null smiles/value", df)

        df = df[df["data_validity_comment"].isna()]
        log.step("2. drop data_validity_comment", df)

        df = df[pd.to_numeric(df[value_col], errors="coerce") > 0]
        log.step("3. drop value <= 0", df)
    else:
        df = df.dropna(subset=["SMILES", pic50_col])
        log.step("1-3. drop null (DDH)", df)

    smi_col = "canonical_smiles" if pic50_col is None else "SMILES"
    std = df[smi_col].apply(standardize)
    df = df.assign(
        clean_smiles=[s[0] for s in std],
        inchikey=[s[1] for s in std],
        n_heavy=[s[2] for s in std],
    )
    df = df.dropna(subset=["clean_smiles", "inchikey"])
    log.step("4-8. RDKit parse + standardize + InChIKey", df)

    if pic50_col is None:
        df["ic50_nM"] = pd.to_numeric(df[value_col])
        df["pIC50"] = 9 - np.log10(df["ic50_nM"])
    else:
        df["pIC50"] = df[pic50_col]
        df["ic50_nM"] = 10 ** (9 - df["pIC50"])
    log.step("9. compute pIC50", df)

    # 10. deduplication by InChIKey
    g = df.groupby("inchikey")
    agg = g.agg(
        clean_smiles=("clean_smiles", "first"),
        molecule_chembl_id=(("molecule_chembl_id" if pic50_col is None
                             else "compound_no"), "first"),
        n_heavy=("n_heavy", "first"),
        pIC50=("pIC50", "median"),
        pIC50_std=("pIC50", "std"),
        n_records=("pIC50", "size"),
    ).reset_index()
    agg["pIC50_std"] = agg["pIC50_std"].fillna(0.0)
    log.step("10a. dedup by InChIKey (median)", agg)

    conflict = agg["pIC50_std"] > DEDUP_SD_MAX
    print(f"      -> conflicting groups (SD > {DEDUP_SD_MAX}): {conflict.sum()}")
    agg = agg[~conflict]
    log.step("10b. drop conflicting groups", agg)

    agg = agg[(agg["pIC50"] >= PIC50_MIN) & (agg["pIC50"] <= PIC50_MAX)]
    log.step(f"11. pIC50 in [{PIC50_MIN}, {PIC50_MAX}]", agg)

    agg = agg[(agg["n_heavy"] >= HEAVY_ATOM_MIN) & (agg["n_heavy"] <= HEAVY_ATOM_MAX)]
    log.step(f"12. heavy atoms in [{HEAVY_ATOM_MIN}, {HEAVY_ATOM_MAX}]", agg)

    agg["ic50_nM"] = 10 ** (9 - agg["pIC50"])
    agg["bioactivity_class"] = agg["pIC50"].apply(bioactivity_class)
    agg["source"] = source
    log.step("13. add bioactivity_class", agg)

    cols = ["molecule_chembl_id", "clean_smiles", "inchikey", "ic50_nM", "pIC50",
            "n_records", "pIC50_std", "n_heavy", "bioactivity_class", "source"]
    return agg[cols].rename(columns={"clean_smiles": "canonical_smiles"}), log.frame()


if __name__ == "__main__":
    print("=== ChEMBL (Mpro) ===")
    chembl_raw = pd.read_csv(DATA_INTERIM / "chembl_after_filter.csv")
    chembl, log_chembl = clean(chembl_raw, source="ChEMBL")

    print("\n=== DDH (COVID-19 Drug Discovery Hackathon) ===")
    ddh_raw = pd.read_csv(DATA_INTERIM / "ddh_labeled.csv")
    ddh, log_ddh = clean(ddh_raw, source="DDH", pic50_col="pIC50")

    # ตรวจ overlap ระหว่างชุดหลักกับ external (ความเสี่ยง R10)
    overlap = set(chembl["inchikey"]) & set(ddh["inchikey"])
    print(f"\n>>> InChIKey overlap ChEMBL vs DDH: {len(overlap)}")
    ddh_ext = ddh[~ddh["inchikey"].isin(overlap)].copy()
    print(f">>> DDH external set after removing overlap: {len(ddh_ext)}")

    chembl.to_csv(DATA_PROC / "clean_bioactivity_master.csv", index=False)
    ddh_ext.to_csv(DATA_PROC / "external_ddh_clean.csv", index=False)
    pd.concat([log_chembl.assign(dataset="ChEMBL"),
               log_ddh.assign(dataset="DDH")]).to_csv(
        DATA_PROC / "cleaning_log.csv", index=False)

    print(f"\nsaved clean_bioactivity_master.csv  shape={chembl.shape}")
    print(chembl["bioactivity_class"].value_counts().to_string())
    print("\npIC50 summary:")
    print(chembl["pIC50"].describe().to_string())

"""Phase 1 — Data Acquisition
ดึง IC50 ของ SARS-CoV-2 Mpro จาก ChEMBL Web Services และเตรียมชุด DDH
"""
import json
import time
import urllib.request
from datetime import datetime

import pandas as pd

from config import (CHEMBL_API, DATA_INTERIM, DATA_RAW, DATASET_ORIG,
                    MPRO_ASSAY_REGEX, TARGET_CHEMBL_ID)

FIELDS = [
    "molecule_chembl_id", "canonical_smiles", "standard_value", "standard_units",
    "standard_relation", "standard_type", "assay_chembl_id", "assay_description",
    "document_chembl_id", "pchembl_value", "data_validity_comment", "activity_id",
]


def _get(url: str, retries: int = 4) -> dict:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception as exc:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1}/{retries} after error: {exc}")
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("unreachable")


def fetch_activities(limit: int = 1000) -> pd.DataFrame:
    url = (f"{CHEMBL_API}/activity?target_chembl_id={TARGET_CHEMBL_ID}"
           f"&standard_type=IC50&standard_units=nM&standard_relation=%3D"
           f"&format=json&limit={limit}")
    rows, page = [], 0
    while url:
        d = _get(url)
        rows.extend(d["activities"])
        page += 1
        total = d["page_meta"]["total_count"]
        print(f"  page {page}: {len(rows)}/{total} records")
        nxt = d["page_meta"]["next"]
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
        time.sleep(0.3)  # เคารพ rate limit ของ EBI
    df = pd.DataFrame(rows)
    return df[[c for c in FIELDS if c in df.columns]]


def prepare_ddh() -> tuple[pd.DataFrame, pd.DataFrame]:
    """อ่าน DDH Data.csv, แยก BLINDED, แปลงหน่วย µM -> pIC50 มาตรฐาน (+6)"""
    p = DATASET_ORIG / "COVID-19 Drug Discovery Data" / "DDH Data.csv"
    df = pd.read_csv(p)
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
    col = "pIC50 (IC50 in microM)"
    df = df.rename(columns={col: "pIC50_uM_basis", "Compound No.": "compound_no"})

    num = pd.to_numeric(df["pIC50_uM_basis"], errors="coerce")
    blinded = df[num.isna()].copy()
    labeled = df[num.notna()].copy()
    labeled["pIC50"] = num[num.notna()] + 6.0  # -log10(uM) -> -log10(M)
    labeled["ic50_nM"] = 10 ** (9 - labeled["pIC50"])
    return labeled, blinded


if __name__ == "__main__":
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    print(f"[{stamp}] Fetching ChEMBL activities for {TARGET_CHEMBL_ID} ...")
    raw = fetch_activities()
    out = DATA_RAW / f"chembl_{TARGET_CHEMBL_ID}_ic50_raw.csv"
    raw.to_csv(out, index=False)
    print(f"saved raw: {out}  shape={raw.shape}")

    mask = raw["assay_description"].str.contains(MPRO_ASSAY_REGEX, case=False,
                                                 na=False, regex=True)
    mpro = raw[mask].copy()
    print(f"assay filter (3CL/Mpro): kept {mask.sum()} / {len(raw)} "
          f"({100 * mask.mean():.1f}%), dropped {(~mask).sum()}")
    mpro.to_csv(DATA_INTERIM / "chembl_after_filter.csv", index=False)

    print("\nตัวอย่าง assay ที่ถูกคัดออก (5 อันดับแรก):")
    for d in raw.loc[~mask, "assay_description"].value_counts().head(5).items():
        print(f"  [{d[1]:>4}] {str(d[0])[:110]}")

    lab, bli = prepare_ddh()
    lab.to_csv(DATA_INTERIM / "ddh_labeled.csv", index=False)
    bli.to_csv(DATA_INTERIM / "ddh_blinded.csv", index=False)
    print(f"\nDDH: labeled={len(lab)}  blinded={len(bli)}")
    print(f"DDH pIC50 (standard) range: {lab['pIC50'].min():.3f} .. "
          f"{lab['pIC50'].max():.3f}  mean={lab['pIC50'].mean():.3f}")

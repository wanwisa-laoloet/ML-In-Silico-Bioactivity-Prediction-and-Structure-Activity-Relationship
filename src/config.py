"""ค่าคงที่กลางของโครงงาน — ทุกสคริปต์ต้อง import จากที่นี่เท่านั้น"""
from pathlib import Path

# ---------- Paths ----------
ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
DATA_INTERIM = ROOT / "data" / "interim"
DATA_PROC = ROOT / "data" / "processed"
DATASET_ORIG = ROOT / "Data_set"
REPORTS = ROOT / "reports"
FIGURES = REPORTS / "figures"
TABLES = REPORTS / "tables"
MODELS = ROOT / "models"

for _p in (DATA_RAW, DATA_INTERIM, DATA_PROC, FIGURES, TABLES, MODELS):
    _p.mkdir(parents=True, exist_ok=True)

# ---------- Reproducibility ----------
RANDOM_STATE = 42

# ---------- ChEMBL ----------
TARGET_CHEMBL_ID = "CHEMBL4523582"  # SARS-CoV-2 Replicase polyprotein 1ab (UniProt P0DTD1)
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
MPRO_ASSAY_REGEX = r"3CL|3C-like|3C like|Mpro|M\s*pro|main\s+protease"

# ---------- Features ----------
FP_RADIUS = 2
FP_BITS = 1024
FP_COLS = [f"FP_{i}" for i in range(FP_BITS)]
DESC_COLS = ["MW", "LogP", "HBD", "HBA"]
FEATURE_COLS = FP_COLS + DESC_COLS  # 1,028 features

# ---------- Cleaning thresholds ----------
PIC50_MIN, PIC50_MAX = 3.0, 11.0
HEAVY_ATOM_MIN, HEAVY_ATOM_MAX = 5, 70
DEDUP_SD_MAX = 1.0  # log units — กลุ่มที่ขัดแย้งกันเกินนี้ให้ตัดทิ้ง

# ---------- Modeling ----------
TEST_SIZE = 0.25
CV_FOLDS = 5
N_JOBS = 4  # ตามสเปกฮาร์ดแวร์ 4 cores


def bioactivity_class(pic50: float) -> str:
    if pic50 >= 6.0:
        return "active"
    if pic50 >= 5.0:
        return "intermediate"
    return "inactive"

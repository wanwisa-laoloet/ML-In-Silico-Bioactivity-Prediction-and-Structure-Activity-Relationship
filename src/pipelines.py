"""Pipeline ที่ป้องกัน Data Leakage — หัวใจของวัตถุประสงค์ O2

StandardScaler ถูก fit ใหม่ทุก train fold ผ่าน ColumnTransformer ภายใน Pipeline
fingerprint bits (0/1) ผ่านแบบ passthrough ไม่ต้อง scale
"""
from collections import defaultdict

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from xgboost import XGBRegressor

from config import DESC_COLS, FP_COLS, N_JOBS, RANDOM_STATE

RDLogger.DisableLog("rdApp.*")


def make_pipe(estimator, fp_cols=FP_COLS, desc_cols=DESC_COLS):
    transformers = []
    if fp_cols:
        transformers.append(("fp", "passthrough", list(fp_cols)))
    if desc_cols:
        transformers.append(("desc", StandardScaler(), list(desc_cols)))
    prep = ColumnTransformer(transformers, remainder="drop")
    return Pipeline([("prep", prep), ("model", estimator)])


def get_estimators():
    """โมเดล 4 ตัวตามข้อกำหนดในหัวข้อ 4.3 ของข้อเสนอโครงงาน"""
    return {
        "Dummy Regressor": DummyRegressor(strategy="mean"),
        "Random Forest": RandomForestRegressor(
            n_estimators=500, random_state=RANDOM_STATE, n_jobs=N_JOBS),
        "XGBoost": XGBRegressor(
            n_estimators=500, learning_rate=0.05, max_depth=6, subsample=0.8,
            colsample_bytree=0.8, random_state=RANDOM_STATE, n_jobs=N_JOBS,
            tree_method="hist", verbosity=0),
        "SVR (RBF)": SVR(kernel="rbf", C=10, gamma="scale", epsilon=0.1),
    }


def scaffold_split(smiles_list, test_frac=0.25):
    """Bemis-Murcko scaffold split — กลุ่มใหญ่เข้า train ก่อน
    โมเลกุลที่มี scaffold เดียวกันจะไม่ถูกแยกไปคนละฝั่ง
    """
    buckets = defaultdict(list)
    for i, s in enumerate(smiles_list):
        try:
            sc = MurckoScaffold.MurckoScaffoldSmiles(smiles=s, includeChirality=False)
        except Exception:  # noqa: BLE001
            sc = ""
        buckets[sc].append(i)
    groups = sorted(buckets.values(), key=len, reverse=True)
    n_train = len(smiles_list) - int(len(smiles_list) * test_frac)
    train, test = [], []
    for g in groups:
        if len(train) + len(g) <= n_train:
            train.extend(g)
        else:
            test.extend(g)
    return np.array(train), np.array(test)

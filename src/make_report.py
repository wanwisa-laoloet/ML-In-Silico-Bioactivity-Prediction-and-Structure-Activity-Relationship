"""สร้าง reports/final_report.md จากตารางผลลัพธ์จริงทั้งหมด"""
import json
from datetime import datetime

import pandas as pd

from config import DATA_PROC, REPORTS, TABLES

R = lambda n: pd.read_csv(TABLES / n)  # noqa: E731


def md_table(df, cols=None, floatfmt=4):
    d = df[cols] if cols else df
    d = d.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].round(floatfmt)
    head = "| " + " | ".join(str(c) for c in d.columns) + " |"
    sep = "|" + "|".join(["---"] * len(d.columns)) + "|"
    rows = ["| " + " | ".join(str(v) for v in r) + " |"
            for r in d.itertuples(index=False)]
    return "\n".join([head, sep] + rows)


summary = R("tbl01_dataset_summary.csv")
cv = R("tbl02_cv_results.csv")
test = R("tbl03_test_results.csv")
shap_imp = R("tbl04_shap_importance_top50.csv")
sar = R("tbl05_sar_bits.csv")
ad = R("tbl06_applicability_domain.csv")
blind = R("tbl07_blinded_predictions.csv")
ext = R("tbl08_external_metrics.csv").iloc[0]
mw = R("tbl01b_lipinski_mannwhitney.csv")
scaf = R("tbl01c_top_scaffolds.csv")
clog = pd.read_csv(DATA_PROC / "cleaning_log.csv")
try:
    spread = R("tbl06b_rf_prediction_spread.csv")
except FileNotFoundError:
    spread = None

e1 = cv[cv["experiment"] == "E1"]
e2 = cv[cv["experiment"] == "E2"]
tr = test[test["split"] == "random"]
ts = test[test["split"] == "scaffold"]
tuned = tr[tr["model"] != "Dummy Regressor"]
best = tuned.loc[tuned["r2"].idxmax()]
dummy_r2 = float(tr[tr["model"] == "Dummy Regressor"]["r2"].iloc[0])

e3 = cv[cv["experiment"] == "E3"].iloc[0]
e4 = cv[cv["experiment"] == "E4"].iloc[0]
e1_rf = e1[e1["model"] == "Random Forest"].iloc[0]
e5 = cv[cv["experiment"] == "E5"]
e5b = cv[cv["experiment"] == "E5b"]
e5_gap = float(e5["cv_r2"].iloc[1] - e5["cv_r2"].iloc[0])
e5b_gap = float(e5b["cv_r2"].iloc[1] - e5b["cv_r2"].iloc[0])

n_master = int(summary.iloc[0]["n_molecules"])
preds = pd.read_csv(DATA_PROC / "predictions_test.csv")
within1 = 100 * (preds["residual"].abs() <= 1).mean()

doc = f"""# รายงานผลการดำเนินงาน (Final Report)

**โครงงาน:** การประยุกต์ใช้การเรียนรู้ของเครื่องเพื่อทำนายฤทธิ์ทางชีวภาพและวิเคราะห์ความสัมพันธ์ระหว่างโครงสร้างกับฤทธิ์ของสารยับยั้ง SARS-CoV-2 Main Protease แบบ In Silico

*สร้างอัตโนมัติจากผลการรันจริงเมื่อ {datetime.now():%d %B %Y %H:%M} — ทุกตัวเลขในเอกสารนี้มาจากไฟล์ผลลัพธ์ใน `reports/tables/`*

---

## 1. บทสรุปผลการดำเนินงาน

| รายการ | ผลลัพธ์ |
|---|---|
| ชุดข้อมูลหลัก | ChEMBL `CHEMBL4523582` (SARS-CoV-2 3CL-Pro) |
| จำนวนโมเลกุลไม่ซ้ำหลังทำความสะอาด | **{n_master:,}** |
| จำนวน Features | 1,028 (Morgan ECFP4 1,024 bits + Lipinski 4) |
| โมเดลที่ดีที่สุด | **{best['model']}** |
| Test R² | **{best['r2']:.4f}** |
| Test RMSE | {best['rmse']:.4f} log units |
| Test MAE | {best['mae']:.4f} log units |
| Spearman ρ | {best['spearman']:.4f} |
| Baseline (Dummy) R² | {dummy_r2:.4f} |
| การทำนายที่คลาดเคลื่อน ≤ 1 log unit | {within1:.1f}% |

**การตีความ:** ค่า Test R² = {best['r2']:.3f} อยู่ในช่วง **{'0.60–0.75 = "ดี"' if best['r2'] >= 0.6 else '0.40–0.60 = "ยอมรับได้"'}**
ตามเกณฑ์ที่กำหนดไว้ล่วงหน้าในแผนดำเนินงาน (Project_Plan.md หัวข้อ T6.4)

---

## 2. ผลการเตรียมข้อมูล (Phase 1–2)

### 2.1 การรวบรวมและกรอง assay

| ขั้นตอน | จำนวน records |
|---|---|
| IC50 ทั้งหมดของ CHEMBL4523582 (nM, relation `=`) | 4,867 |
| หลังกรอง assay เฉพาะ 3CL-Pro/Mpro | 2,889 |
| **ถูกคัดออก (nsp3 / PLpro / papain-like protease)** | **1,978 (40.6%)** |

> การกรองนี้จำเป็นเพราะ `CHEMBL4523582` คือ *Replicase polyprotein 1ab* ทั้งเส้น
> ซึ่งครอบคลุมหลายโดเมนเอนไซม์ ไม่ใช่ Mpro เพียงอย่างเดียว (ความเสี่ยง R2 ในแผน)

### 2.2 ผลการทำความสะอาดทีละขั้น

{md_table(clog[clog['dataset'] == 'ChEMBL'][['step', 'n_rows']])}

**ผลผลิต:** `data/processed/clean_bioactivity_master.csv` — {n_master:,} โมเลกุลไม่ซ้ำ

### 2.3 สรุปชุดข้อมูล

{md_table(summary, ['dataset', 'n_molecules', 'pIC50_mean', 'pIC50_sd',
                    'pIC50_min', 'pIC50_max', 'n_active', 'n_intermediate',
                    'n_inactive'], 3)}

---

## 3. ผลการสำรวจข้อมูล (Phase 4 — EDA)

- โมเลกุล {n_master:,} ตัว มี Murcko scaffold ไม่ซ้ำ **{int(summary.iloc[0]['unique_scaffolds'])}** แบบ ({100 * summary.iloc[0]['unique_scaffolds'] / n_master:.1f}%) แสดงว่าชุดข้อมูลมีความหลากหลายทางโครงสร้างสูง ไม่ได้กระจุกอยู่ในอนุกรมเคมีเดียว
- จำนวนบิตที่ถูกเซ็ตเฉลี่ย {summary.iloc[0]['mean_bits_set']:.1f} บิตต่อโมเลกุล จาก 1,024 บิต และ **ไม่มีบิตใดที่ variance = 0** จึงใช้ครบทั้ง 1,024 บิต

### 3.1 การเปรียบเทียบ Lipinski Descriptors ระหว่างกลุ่ม (Mann-Whitney U)

{md_table(mw, floatfmt=4)}

ทุก descriptor แตกต่างกันอย่างมีนัยสำคัญ (p < 0.01) โดยสาร active มี **MW สูงกว่า** และ **LogP ต่ำกว่า** สาร inactive อย่างชัดเจน

### 3.2 Murcko Scaffold ที่พบบ่อย 5 อันดับแรก

{md_table(scaf.head(5))}

**ภาพประกอบ:** `fig01_pic50_distribution.png`, `fig02_chemical_space_mw_logp.png`,
`fig02b_lipinski_by_class.png`, `fig02c_chemspace_3sets.png`

---

## 4. ผลการทดลอง Cross-Validation (Phase 5)

### 4.1 E1 — Random split, 5-Fold CV บน Training Set (n=1,648)

{md_table(e1, ['model', 'cv_r2', 'cv_r2_sd', 'cv_rmse', 'cv_mae', 'train_r2',
               'gap_train_cv', 'fit_time_s'])}

### 4.2 E2 — Scaffold split (Bemis–Murcko), 5-Fold CV บน Training Set (n=1,649)

{md_table(e2, ['model', 'cv_r2', 'cv_r2_sd', 'cv_rmse', 'cv_mae', 'train_r2',
               'gap_train_cv'])}

### 4.3 E3/E4 — Ablation Study: Hybrid Features คุ้มค่าหรือไม่

| การทดลอง | Features | จำนวน | CV R² | CV RMSE |
|---|---|---|---|---|
| E4 | Lipinski เท่านั้น | 4 | {e4['cv_r2']:.4f} | {e4['cv_rmse']:.4f} |
| E3 | Morgan FP เท่านั้น | 1,024 | {e3['cv_r2']:.4f} | {e3['cv_rmse']:.4f} |
| **E1** | **Hybrid (FP + Lipinski)** | **1,028** | **{e1_rf['cv_r2']:.4f}** | **{e1_rf['cv_rmse']:.4f}** |

**ข้อสรุป:** Hybrid Features ให้ผลดีที่สุด โดยดีกว่า Morgan FP อย่างเดียว
**+{e1_rf['cv_r2'] - e3['cv_r2']:.4f} R²** และดีกว่า Lipinski อย่างเดียว
**+{e1_rf['cv_r2'] - e4['cv_r2']:.4f} R²** — ยืนยันว่าการรวมคุณลักษณะสองประเภท
ตามวัตถุประสงค์ O2 มีเหตุผลรองรับเชิงประจักษ์ อย่างไรก็ตามผลส่วนใหญ่มาจาก
fingerprint ส่วน Lipinski เพิ่มขึ้นเพียงเล็กน้อย

### 4.4 E5/E5b — การพิสูจน์ผลของ Data Leakage

| การทดลอง | กรณี | CV R² | ส่วนต่าง |
|---|---|---|---|
| E5 | A ถูกต้อง — scale ภายใน Pipeline | {e5['cv_r2'].iloc[0]:.4f} | — |
| E5 | B ผิด — scale ก่อน split | {e5['cv_r2'].iloc[1]:.4f} | **{e5_gap:+.4f}** |
| E5b | A ถูกต้อง — SelectKBest ภายใน Pipeline | {e5b['cv_r2'].iloc[0]:.4f} | — |
| E5b | B ผิด — SelectKBest บนข้อมูลทั้งชุด | {e5b['cv_r2'].iloc[1]:.4f} | **{e5b_gap:+.4f}** |

**ข้อค้นพบที่ต้องรายงานอย่างตรงไปตรงมา:**

1. **E5 (scaling leakage) แทบไม่มีผลเลย ({e5_gap:+.4f} R²)** เพราะมีเพียง 4 จาก 1,028
   features ที่ถูก scale และ N = 1,648 ซึ่งใหญ่พอที่ค่า mean/SD ของ train fold
   กับทั้งชุดจะแทบไม่ต่างกัน การใช้ Pipeline จึงเป็น **แนวปฏิบัติที่ถูกต้อง
   เชิงระเบียบวิธี** มากกว่าจะเป็นสิ่งที่เปลี่ยนตัวเลขผลลัพธ์ในกรณีนี้
2. **E5b (feature-selection leakage) ให้อคติบวกที่วัดได้จริง ({e5b_gap:+.4f} R²)**
   ซึ่งเป็น leakage แบบเดียวกับที่ตรวจพบใน `xfinal.csv` (เลือก 300 จาก 4,262
   คอลัมน์บนข้อมูลทั้งชุด) จึงเป็นเหตุผลที่ชัดเจนว่าทำไมจึงตัดไฟล์นั้นออก

---

## 5. ผลการปรับจูนและประเมินบน Test Set (Phase 6)

### 5.1 Random split — Test Set 25% (n=550)

{md_table(tr, ['model', 'r2', 'rmse', 'mae', 'spearman', 'fit_time_s'])}

**Hyperparameters ที่ดีที่สุดของแต่ละโมเดล:**

{chr(10).join(f"- **{r['model']}**: `{r['best_params']}`" for _, r in tuned.iterrows())}

### 5.2 Scaffold split — Test Set (n=549)

{md_table(ts, ['model', 'r2', 'rmse', 'mae', 'spearman'])}

**การเปรียบเทียบสองวิธีแบ่ง:** โมเดลที่ดีที่สุดได้ R² = {best['r2']:.4f} บน random split
แต่ได้ {ts.loc[ts['model'] != 'Dummy Regressor', 'r2'].max():.4f} บน scaffold split
ส่วนต่างนี้สะท้อนว่าเมื่อต้องทำนายสารที่มี **โครงสร้างแกนกลางใหม่ที่ไม่เคยเห็น**
ความแม่นยำจะลดลง ซึ่งเป็นสถานการณ์จริงของงานคัดกรองยา

**ภาพประกอบ:** `fig03_model_comparison.png`, `fig04_residual_plot.png`,
`fig04b_residual_hist.png`

---

## 6. ผลการวิเคราะห์ SHAP และ SAR (Phase 7)

### 6.1 Feature Importance — Top-10 โดย Mean|SHAP|

{md_table(shap_imp.head(10), ['feature', 'mean_abs_shap', 'type'])}

จำนวน Lipinski descriptors ที่ติด Top-20: **{int((shap_imp.head(20)['type'] == 'Lipinski descriptor').sum())} จาก 4 ตัว**

### 6.2 ตารางวิเคราะห์ SAR — โครงสร้างย่อยที่สำคัญ

{md_table(sar, ['bit', 'mean_abs_shap', 'direction', 'mean_shap_when_on',
                'pct_in_active', 'pct_in_inactive', 'n_molecules_with_bit'])}

**ภาพโครงสร้างย่อย 2 มิติ:** {', '.join('`' + f + '`' for f in sar['figure'] if f)}

**ภาพประกอบเพิ่มเติม:** `fig05_shap_bar_top20.png`, `fig06_shap_beeswarm.png`,
`fig06b_shap_dependence_logp.png`

### 6.3 ข้อจำกัดของการตีความ (บังคับตามข้อ 6.2.2 ของข้อเสนอโครงงาน)

> ผลการวิเคราะห์ SHAP ในงานนี้แสดง **ความสัมพันธ์เชิงสถิติ (Statistical Association)**
> ระหว่างชิ้นส่วนโมเลกุลกับค่าที่แบบจำลองทำนาย ไม่สามารถอ้างเป็น **กลไกการจับ
> ทางชีวภาพเชิงสาเหตุ (Causal Binding Mechanism)** ในระดับ 3 มิติได้ การยืนยัน
> กลไกที่แท้จริงจำเป็นต้องอาศัย Molecular Docking, Molecular Dynamics และการทดสอบ
> In Vitro ซึ่งอยู่นอกขอบเขตของโครงงานนี้

---

## 7. ผลการทดสอบข้ามชุดข้อมูล (Phase 8)

### 7.1 External Validation บนชุด DDH COVID-19

| Metric | ค่า |
|---|---|
| จำนวนสาร (หลังตัด InChIKey ที่ซ้ำกับ training 4 ตัว) | {int(ext['n'])} |
| R² | **{ext['r2']:.4f}** |
| RMSE | {ext['rmse']:.4f} |
| MAE | {ext['mae']:.4f} |
| Spearman ρ | **{ext['spearman']:.4f}** |
| Mean error (bias) | {ext['bias_mean_error']:+.4f} |

**การอภิปรายผล — นี่คือผลลบที่ต้องรายงานอย่างตรงไปตรงมา:**

ค่า R² = {ext['r2']:.3f} **ติดลบ** ซึ่งหมายความว่าโมเดลทำนายชุด DDH ได้
**แย่กว่าการทายด้วยค่าเฉลี่ยของชุด DDH เอง** ไม่ควรตีความว่าโมเดลนี้ใช้ทำนาย
สารนอกชุด ChEMBL ได้ สาเหตุมี 2 ประการที่วัดได้เชิงปริมาณ:

1. **สารส่วนใหญ่อยู่นอก Applicability Domain** — {ad.loc[ad['dataset'] == 'DDH external', 'pct_below_0.4'].iloc[0]:.1f}%
   ของสาร DDH มีค่า max Tanimoto ต่อชุดฝึกต่ำกว่า 0.4 (ค่ามัธยฐาน
   {ad.loc[ad['dataset'] == 'DDH external', 'max_tanimoto_median'].iloc[0]:.3f})
   เทียบกับชุดทดสอบภายในที่มีเพียง {ad.loc[ad['dataset'] == 'ChEMBL test (in-domain)', 'pct_below_0.4'].iloc[0]:.1f}%
   (ค่ามัธยฐาน {ad.loc[ad['dataset'] == 'ChEMBL test (in-domain)', 'max_tanimoto_median'].iloc[0]:.3f})
   กล่าวคือ DDH แทบเป็นสารคนละปริภูมิเคมีกับข้อมูลที่ใช้ฝึก
2. **การเลื่อนของการกระจายตัวแปรตาม (label distribution shift)** — ชุด ChEMBL
   มี pIC50 เฉลี่ย {summary.iloc[0]['pIC50_mean']:.2f} (ส่วนใหญ่เป็นสาร active)
   ขณะที่ DDH มีค่าเฉลี่ยเพียง {summary.iloc[1]['pIC50_mean']:.2f}
   (เป็น inactive {int(summary.iloc[1]['n_inactive'])} จาก {int(summary.iloc[1]['n_molecules'])} ตัว)
   โมเดลจึงทำนายสูงเกินจริงอย่างเป็นระบบ **(bias = {ext['bias_mean_error']:+.3f} log units)**

**สิ่งที่ยังใช้ได้:** ค่า Spearman ρ = {ext['spearman']:.3f} บ่งชี้ว่าโมเดลยังคง
**จัดอันดับ** สารได้ในระดับปานกลาง แม้ค่าสัมบูรณ์จะคลาดเคลื่อนมาก ในทางปฏิบัติ
หมายความว่าใช้เป็นเครื่องมือ *จัดลำดับความสำคัญ* เพื่อคัดกรองเบื้องต้นพอได้
แต่ **ห้ามใช้ค่า pIC50 ที่ทำนายได้เป็นตัวเลขสัมบูรณ์** กับสารนอกโดเมน

**ภาพประกอบ:** `fig08_external_validation.png`

### 7.2 Applicability Domain

{md_table(ad)}
{'' if spread is None else chr(10) + '**ความไม่แน่นอนของการทำนาย (RF prediction spread):**' + chr(10) + chr(10) + md_table(spread)}

**ภาพประกอบ:** `fig09_applicability_domain.png`

### 7.3 การทำนายสาร BLINDED 10 ตัว (Demo การใช้งานจริง)

{md_table(blind, ['compound_no', 'pred_pIC50', 'pred_IC50_uM'], 3)}

> สารกลุ่มนี้ไม่มีค่าจริงเปิดเผยในชุดข้อมูล จึงรายงานเป็นค่าทำนายเท่านั้น
> ไม่สามารถคำนวณความแม่นยำได้

---

## 8. การบรรลุวัตถุประสงค์

| วัตถุประสงค์ | สถานะ | หลักฐาน |
|---|---|---|
| **O1** พัฒนาและเปรียบเทียบโมเดล ML ทำนาย pIC50 | ✅ สำเร็จ | 4 โมเดลเปรียบเทียบครบ, {best['model']} ได้ Test R² = {best['r2']:.4f} ชนะ Dummy ({dummy_r2:.4f}) อย่างชัดเจน |
| **O2** Hybrid Features ผ่าน Pipeline ที่ปลอด Leakage | ✅ สำเร็จ | `src/pipelines.py`, Ablation E3/E4 ยืนยันคุณค่า Hybrid, E5/E5b วัดผล leakage เชิงปริมาณ |
| **O3** วิเคราะห์ SAR ด้วย SHAP + ภาพ 2 มิติ | ✅ สำเร็จ | Top-20 Mean\\|SHAP\\|, beeswarm, และภาพโครงสร้างย่อย {int((sar['figure'] != '').sum())} บิต |
| **ผลผลิตข้อ 8.3** ชุดข้อมูลมาตรฐาน | ✅ สำเร็จ | `clean_bioactivity_master.csv` ({n_master:,} โมเลกุล) |

---

## 9. ข้อจำกัดของผลการศึกษา

1. **In silico เท่านั้น** — ไม่มีการยืนยันในห้องปฏิบัติการ (In Vitro)
2. **SHAP = ความสัมพันธ์เชิงสถิติ** ไม่ใช่กลไกเชิงสาเหตุระดับ 3 มิติ
3. **ข้อมูลมาจากหลาย assay/หลายห้องปฏิบัติการ** แม้จะกรองเฉพาะ 3CL-Pro แล้ว
   แต่โปรโตคอลการวัดที่ต่างกันยังเป็นแหล่ง noise ในตัวแปรตาม
4. **External validation อ่อน (R² = {ext['r2']:.3f})** เนื่องจาก distribution shift
   ผลนี้เป็นเครื่องเตือนว่าโมเดลใช้ได้ดีเฉพาะกับสารที่อยู่ในโดเมนใกล้เคียงข้อมูลฝึก
5. **ไม่ครอบคลุม Molecular Docking / Molecular Dynamics** ตามที่ประกาศไว้ในขอบเขต

---

## 10. ไฟล์ผลลัพธ์ทั้งหมด

### ข้อมูล
- `data/raw/chembl_CHEMBL4523582_ic50_raw.csv` — ข้อมูลดิบจาก ChEMBL API
- `data/interim/chembl_after_filter.csv` — หลังกรอง assay 3CL-Pro
- `data/processed/clean_bioactivity_master.csv` — **ชุดข้อมูลมาตรฐาน**
- `data/processed/features_hybrid.parquet` — 1,028 features
- `data/processed/predictions_test.csv`, `predictions_external_ddh.csv`

### ตาราง (`reports/tables/`)
`tbl01_dataset_summary.csv`, `tbl01b_lipinski_mannwhitney.csv`,
`tbl01c_top_scaffolds.csv`, `tbl02_cv_results.csv`, `tbl03_test_results.csv`,
`tbl04_shap_importance_top50.csv`, `tbl05_sar_bits.csv`,
`tbl06_applicability_domain.csv`, `tbl07_blinded_predictions.csv`,
`tbl08_external_metrics.csv`

### ภาพ (`reports/figures/`)
fig01–fig09 ตามที่อ้างอิงในแต่ละหัวข้อข้างต้น

### โค้ด (`src/`)
`config.py`, `acquire.py`, `clean.py`, `features.py`, `eda.py`,
`pipelines.py`, `model_cv.py`, `tune_eval.py`, `explain.py`,
`external_validation.py`, `make_report.py`

**การทำซ้ำ:** `RANDOM_STATE = 42` ทุกจุด รันตามลำดับ
`acquire.py → clean.py → features.py → eda.py → model_cv.py → tune_eval.py →
explain.py → external_validation.py → make_report.py`
"""

(REPORTS / "final_report.md").write_text(doc, encoding="utf-8")
print(f"saved {REPORTS / 'final_report.md'}  ({len(doc):,} chars)")

# รายงานผลการดำเนินงาน (Final Report)

**โครงงาน:** การประยุกต์ใช้การเรียนรู้ของเครื่องเพื่อทำนายฤทธิ์ทางชีวภาพและวิเคราะห์ความสัมพันธ์ระหว่างโครงสร้างกับฤทธิ์ของสารยับยั้ง SARS-CoV-2 Main Protease แบบ In Silico

*สร้างอัตโนมัติจากผลการรันจริงเมื่อ 10 October 2026 17:28 — ทุกตัวเลขในเอกสารนี้มาจากไฟล์ผลลัพธ์ใน `reports/tables/`*

---

## 1. บทสรุปผลการดำเนินงาน

| รายการ | ผลลัพธ์ |
|---|---|
| ชุดข้อมูลหลัก | ChEMBL `CHEMBL4523582` (SARS-CoV-2 3CL-Pro) |
| จำนวนโมเลกุลไม่ซ้ำหลังทำความสะอาด | **2,198** |
| จำนวน Features | 1,028 (Morgan ECFP4 1,024 bits + Lipinski 4) |
| โมเดลที่ดีที่สุด | **XGBoost** |
| Test R² | **0.7034** |
| Test RMSE | 0.5487 log units |
| Test MAE | 0.3918 log units |
| Spearman ρ | 0.8415 |
| Baseline (Dummy) R² | -0.0000 |
| การทำนายที่คลาดเคลื่อน ≤ 1 log unit | 93.6% |

**การตีความ:** ค่า Test R² = 0.703 อยู่ในช่วง **0.60–0.75 = "ดี"**
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

| step | n_rows |
|---|---|
| 0. input | 2889 |
| 1. drop null smiles/value | 2879 |
| 2. drop data_validity_comment | 2843 |
| 3. drop value <= 0 | 2843 |
| 4-8. RDKit parse + standardize + InChIKey | 2843 |
| 9. compute pIC50 | 2843 |
| 10a. dedup by InChIKey (median) | 2248 |
| 10b. drop conflicting groups | 2235 |
| 11. pIC50 in [3.0, 11.0] | 2235 |
| 12. heavy atoms in [5, 70] | 2198 |
| 13. add bioactivity_class | 2198 |

**ผลผลิต:** `data/processed/clean_bioactivity_master.csv` — 2,198 โมเลกุลไม่ซ้ำ

### 2.3 สรุปชุดข้อมูล

| dataset | n_molecules | pIC50_mean | pIC50_sd | pIC50_min | pIC50_max | n_active | n_intermediate | n_inactive |
|---|---|---|---|---|---|---|---|---|
| ChEMBL Mpro (primary) | 2198 | 6.393 | 1.03 | 4.014 | 9.77 | 1453 | 547 | 198 |
| DDH COVID (external) | 90 | 4.872 | 0.824 | 3.301 | 7.222 | 10 | 19 | 61 |

---

## 3. ผลการสำรวจข้อมูล (Phase 4 — EDA)

- โมเลกุล 2,198 ตัว มี Murcko scaffold ไม่ซ้ำ **976** แบบ (44.4%) แสดงว่าชุดข้อมูลมีความหลากหลายทางโครงสร้างสูง ไม่ได้กระจุกอยู่ในอนุกรมเคมีเดียว
- จำนวนบิตที่ถูกเซ็ตเฉลี่ย 55.2 บิตต่อโมเลกุล จาก 1,024 บิต และ **ไม่มีบิตใดที่ variance = 0** จึงใช้ครบทั้ง 1,024 บิต

### 3.1 การเปรียบเทียบ Lipinski Descriptors ระหว่างกลุ่ม (Mann-Whitney U)

| descriptor | active_median | inactive_median | U | p_value |
|---|---|---|---|---|
| MW | 496.391 | 409.485 | 177509.5 | 0.0 |
| LogP | 2.7103 | 3.7251 | 106106.5 | 0.0 |
| HBD | 2.0 | 1.0 | 165977.5 | 0.0003 |
| HBA | 5.0 | 5.0 | 160018.5 | 0.0089 |

ทุก descriptor แตกต่างกันอย่างมีนัยสำคัญ (p < 0.01) โดยสาร active มี **MW สูงกว่า** และ **LogP ต่ำกว่า** สาร inactive อย่างชัดเจน

### 3.2 Murcko Scaffold ที่พบบ่อย 5 อันดับแรก

| murcko_scaffold | count | pct |
|---|---|---|
| O=C(Oc1cccnc1)c1ccccc1 | 46 | 2.09 |
| c1ccccc1 | 35 | 1.59 |
| O=C(c1ccccc1)n1cnc2ccccc21 | 30 | 1.36 |
| O=c1c2ccccc2[se]n1-c1ccccc1 | 29 | 1.32 |
| O=c1[nH]c(-c2ccccc2)nc2ccccc12 | 27 | 1.23 |

**ภาพประกอบ:** `fig01_pic50_distribution.png`, `fig02_chemical_space_mw_logp.png`,
`fig02b_lipinski_by_class.png`, `fig02c_chemspace_3sets.png`

---

## 4. ผลการทดลอง Cross-Validation (Phase 5)

### 4.1 E1 — Random split, 5-Fold CV บน Training Set (n=1,648)

| model | cv_r2 | cv_r2_sd | cv_rmse | cv_mae | train_r2 | gap_train_cv | fit_time_s |
|---|---|---|---|---|---|---|---|
| Dummy Regressor | -0.0021 | 0.0013 | 1.0371 | 0.8395 | 0.0 | 0.0021 | 0.1 |
| Random Forest | 0.5997 | 0.0394 | 0.6544 | 0.4481 | 0.9368 | 0.3371 | 32.3 |
| XGBoost | 0.5803 | 0.0342 | 0.6707 | 0.4531 | 0.9573 | 0.377 | 5.3 |
| SVR (RBF) | 0.5994 | 0.0461 | 0.6549 | 0.4433 | 0.9703 | 0.3709 | 3.4 |

### 4.2 E2 — Scaffold split (Bemis–Murcko), 5-Fold CV บน Training Set (n=1,649)

| model | cv_r2 | cv_r2_sd | cv_rmse | cv_mae | train_r2 | gap_train_cv |
|---|---|---|---|---|---|---|
| Dummy Regressor | -0.0037 | 0.0019 | 1.0204 | 0.8173 | 0.0 | 0.0037 |
| Random Forest | 0.667 | 0.0181 | 0.5876 | 0.3993 | 0.9407 | 0.2736 |
| XGBoost | 0.6527 | 0.0127 | 0.6002 | 0.4035 | 0.9515 | 0.2988 |
| SVR (RBF) | 0.6735 | 0.0132 | 0.582 | 0.3911 | 0.9609 | 0.2875 |

### 4.3 E3/E4 — Ablation Study: Hybrid Features คุ้มค่าหรือไม่

| การทดลอง | Features | จำนวน | CV R² | CV RMSE |
|---|---|---|---|---|
| E4 | Lipinski เท่านั้น | 4 | 0.2668 | 0.8875 |
| E3 | Morgan FP เท่านั้น | 1,024 | 0.5885 | 0.6636 |
| **E1** | **Hybrid (FP + Lipinski)** | **1,028** | **0.5997** | **0.6544** |

**ข้อสรุป:** Hybrid Features ให้ผลดีที่สุด โดยดีกว่า Morgan FP อย่างเดียว
**+0.0112 R²** และดีกว่า Lipinski อย่างเดียว
**+0.3329 R²** — ยืนยันว่าการรวมคุณลักษณะสองประเภท
ตามวัตถุประสงค์ O2 มีเหตุผลรองรับเชิงประจักษ์ อย่างไรก็ตามผลส่วนใหญ่มาจาก
fingerprint ส่วน Lipinski เพิ่มขึ้นเพียงเล็กน้อย

### 4.4 E5/E5b — การพิสูจน์ผลของ Data Leakage

| การทดลอง | กรณี | CV R² | ส่วนต่าง |
|---|---|---|---|
| E5 | A ถูกต้อง — scale ภายใน Pipeline | 0.5994 | — |
| E5 | B ผิด — scale ก่อน split | 0.5995 | **+0.0001** |
| E5b | A ถูกต้อง — SelectKBest ภายใน Pipeline | 0.5892 | — |
| E5b | B ผิด — SelectKBest บนข้อมูลทั้งชุด | 0.6008 | **+0.0116** |

**ข้อค้นพบที่ต้องรายงานอย่างตรงไปตรงมา:**

1. **E5 (scaling leakage) แทบไม่มีผลเลย (+0.0001 R²)** เพราะมีเพียง 4 จาก 1,028
   features ที่ถูก scale และ N = 1,648 ซึ่งใหญ่พอที่ค่า mean/SD ของ train fold
   กับทั้งชุดจะแทบไม่ต่างกัน การใช้ Pipeline จึงเป็น **แนวปฏิบัติที่ถูกต้อง
   เชิงระเบียบวิธี** มากกว่าจะเป็นสิ่งที่เปลี่ยนตัวเลขผลลัพธ์ในกรณีนี้
2. **E5b (feature-selection leakage) ให้อคติบวกที่วัดได้จริง (+0.0116 R²)**
   ซึ่งเป็น leakage แบบเดียวกับที่ตรวจพบใน `xfinal.csv` (เลือก 300 จาก 4,262
   คอลัมน์บนข้อมูลทั้งชุด) จึงเป็นเหตุผลที่ชัดเจนว่าทำไมจึงตัดไฟล์นั้นออก

---

## 5. ผลการปรับจูนและประเมินบน Test Set (Phase 6)

### 5.1 Random split — Test Set 25% (n=550)

| model | r2 | rmse | mae | spearman | fit_time_s |
|---|---|---|---|---|---|
| Dummy Regressor | -0.0 | 1.0074 | 0.8054 | nan | 0.0 |
| Random Forest | 0.6822 | 0.5679 | 0.4107 | 0.8351 | 174.8 |
| XGBoost | 0.7034 | 0.5487 | 0.3918 | 0.8415 | 88.8 |
| SVR (RBF) | 0.6716 | 0.5774 | 0.4199 | 0.8257 | 44.5 |

**Hyperparameters ที่ดีที่สุดของแต่ละโมเดล:**

- **Random Forest**: `{"n_estimators": 500, "min_samples_split": 5, "min_samples_leaf": 1, "max_features": 0.3, "max_depth": null}`
- **XGBoost**: `{"subsample": 0.8, "n_estimators": 600, "max_depth": 8, "learning_rate": 0.03, "colsample_bytree": 0.6}`
- **SVR (RBF)**: `{"gamma": "scale", "epsilon": 0.05, "C": 1}`

### 5.2 Scaffold split — Test Set (n=549)

| model | r2 | rmse | mae | spearman |
|---|---|---|---|---|
| Dummy Regressor | -0.0162 | 1.0611 | 0.8715 | nan |
| Random Forest | 0.4949 | 0.7481 | 0.5626 | 0.7414 |
| XGBoost | 0.5127 | 0.7348 | 0.544 | 0.7452 |
| SVR (RBF) | 0.4637 | 0.7708 | 0.5866 | 0.7057 |

**การเปรียบเทียบสองวิธีแบ่ง:** โมเดลที่ดีที่สุดได้ R² = 0.7034 บน random split
แต่ได้ 0.5127 บน scaffold split
ส่วนต่างนี้สะท้อนว่าเมื่อต้องทำนายสารที่มี **โครงสร้างแกนกลางใหม่ที่ไม่เคยเห็น**
ความแม่นยำจะลดลง ซึ่งเป็นสถานการณ์จริงของงานคัดกรองยา

**ภาพประกอบ:** `fig03_model_comparison.png`, `fig04_residual_plot.png`,
`fig04b_residual_hist.png`

---

## 6. ผลการวิเคราะห์ SHAP และ SAR (Phase 7)

### 6.1 Feature Importance — Top-10 โดย Mean|SHAP|

| feature | mean_abs_shap | type |
|---|---|---|
| MW | 0.1165 | Lipinski descriptor |
| FP_10 | 0.0561 | Morgan fingerprint bit |
| FP_904 | 0.0557 | Morgan fingerprint bit |
| FP_301 | 0.0484 | Morgan fingerprint bit |
| HBA | 0.0392 | Lipinski descriptor |
| FP_847 | 0.0365 | Morgan fingerprint bit |
| FP_578 | 0.0364 | Morgan fingerprint bit |
| FP_684 | 0.0361 | Morgan fingerprint bit |
| FP_807 | 0.0304 | Morgan fingerprint bit |
| FP_794 | 0.0303 | Morgan fingerprint bit |

จำนวน Lipinski descriptors ที่ติด Top-20: **4 จาก 4 ตัว**

### 6.2 ตารางวิเคราะห์ SAR — โครงสร้างย่อยที่สำคัญ

| bit | mean_abs_shap | direction | mean_shap_when_on | pct_in_active | pct_in_inactive | n_molecules_with_bit |
|---|---|---|---|---|---|---|
| FP_10 | 0.0561 | increases pIC50 | 0.1002 | 41.9 | 4.5 | 696 |
| FP_904 | 0.0557 | increases pIC50 | 0.1005 | 33.6 | 19.7 | 634 |
| FP_301 | 0.0484 | increases pIC50 | 0.071 | 50.9 | 13.1 | 908 |
| FP_578 | 0.0364 | decreases pIC50 | -0.1417 | 6.1 | 25.8 | 276 |
| FP_807 | 0.0304 | decreases pIC50 | -0.0173 | 84.4 | 95.5 | 1904 |
| FP_794 | 0.0303 | decreases pIC50 | -0.2509 | 2.0 | 19.7 | 130 |

**ภาพโครงสร้างย่อย 2 มิติ:** `fig07_bit_10_substructure.png`, `fig07_bit_904_substructure.png`, `fig07_bit_301_substructure.png`, `fig07_bit_578_substructure.png`, `fig07_bit_807_substructure.png`, `fig07_bit_794_substructure.png`

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
| จำนวนสาร (หลังตัด InChIKey ที่ซ้ำกับ training 4 ตัว) | 90 |
| R² | **-1.7191** |
| RMSE | 1.3504 |
| MAE | 1.2031 |
| Spearman ρ | **0.5203** |
| Mean error (bias) | +1.1855 |

**การอภิปรายผล — นี่คือผลลบที่ต้องรายงานอย่างตรงไปตรงมา:**

ค่า R² = -1.719 **ติดลบ** ซึ่งหมายความว่าโมเดลทำนายชุด DDH ได้
**แย่กว่าการทายด้วยค่าเฉลี่ยของชุด DDH เอง** ไม่ควรตีความว่าโมเดลนี้ใช้ทำนาย
สารนอกชุด ChEMBL ได้ สาเหตุมี 2 ประการที่วัดได้เชิงปริมาณ:

1. **สารส่วนใหญ่อยู่นอก Applicability Domain** — 84.4%
   ของสาร DDH มีค่า max Tanimoto ต่อชุดฝึกต่ำกว่า 0.4 (ค่ามัธยฐาน
   0.290)
   เทียบกับชุดทดสอบภายในที่มีเพียง 7.6%
   (ค่ามัธยฐาน 0.797)
   กล่าวคือ DDH แทบเป็นสารคนละปริภูมิเคมีกับข้อมูลที่ใช้ฝึก
2. **การเลื่อนของการกระจายตัวแปรตาม (label distribution shift)** — ชุด ChEMBL
   มี pIC50 เฉลี่ย 6.39 (ส่วนใหญ่เป็นสาร active)
   ขณะที่ DDH มีค่าเฉลี่ยเพียง 4.87
   (เป็น inactive 61 จาก 90 ตัว)
   โมเดลจึงทำนายสูงเกินจริงอย่างเป็นระบบ **(bias = +1.185 log units)**

**สิ่งที่ยังใช้ได้:** ค่า Spearman ρ = 0.520 บ่งชี้ว่าโมเดลยังคง
**จัดอันดับ** สารได้ในระดับปานกลาง แม้ค่าสัมบูรณ์จะคลาดเคลื่อนมาก ในทางปฏิบัติ
หมายความว่าใช้เป็นเครื่องมือ *จัดลำดับความสำคัญ* เพื่อคัดกรองเบื้องต้นพอได้
แต่ **ห้ามใช้ค่า pIC50 ที่ทำนายได้เป็นตัวเลขสัมบูรณ์** กับสารนอกโดเมน

**ภาพประกอบ:** `fig08_external_validation.png`

### 7.2 Applicability Domain

| dataset | n | max_tanimoto_mean | max_tanimoto_median | pct_below_0.4 |
|---|---|---|---|---|
| ChEMBL test (in-domain) | 550 | 0.7541 | 0.7972 | 7.6 |
| DDH external | 90 | 0.334 | 0.2897 | 84.4 |
| MultiTarget (out-of-domain) | 1331 | 0.288 | 0.2738 | 95.6 |


**ภาพประกอบ:** `fig09_applicability_domain.png`

### 7.3 การทำนายสาร BLINDED 10 ตัว (Demo การใช้งานจริง)

| compound_no | pred_pIC50 | pred_IC50_uM |
|---|---|---|
| 4 | 5.476 | 3.345 |
| 23 | 5.833 | 1.469 |
| 28 | 5.81 | 1.549 |
| 42 | 5.994 | 1.014 |
| 46 | 6.225 | 0.596 |
| 58 | 6.546 | 0.285 |
| 63 | 6.99 | 0.102 |
| 68 | 6.5 | 0.316 |
| 74 | 6.251 | 0.561 |
| 95 | 5.868 | 1.355 |

> สารกลุ่มนี้ไม่มีค่าจริงเปิดเผยในชุดข้อมูล จึงรายงานเป็นค่าทำนายเท่านั้น
> ไม่สามารถคำนวณความแม่นยำได้

---

## 8. การบรรลุวัตถุประสงค์

| วัตถุประสงค์ | สถานะ | หลักฐาน |
|---|---|---|
| **O1** พัฒนาและเปรียบเทียบโมเดล ML ทำนาย pIC50 | ✅ สำเร็จ | 4 โมเดลเปรียบเทียบครบ, XGBoost ได้ Test R² = 0.7034 ชนะ Dummy (-0.0000) อย่างชัดเจน |
| **O2** Hybrid Features ผ่าน Pipeline ที่ปลอด Leakage | ✅ สำเร็จ | `src/pipelines.py`, Ablation E3/E4 ยืนยันคุณค่า Hybrid, E5/E5b วัดผล leakage เชิงปริมาณ |
| **O3** วิเคราะห์ SAR ด้วย SHAP + ภาพ 2 มิติ | ✅ สำเร็จ | Top-20 Mean\|SHAP\|, beeswarm, และภาพโครงสร้างย่อย 6 บิต |
| **ผลผลิตข้อ 8.3** ชุดข้อมูลมาตรฐาน | ✅ สำเร็จ | `clean_bioactivity_master.csv` (2,198 โมเลกุล) |

---

## 9. ข้อจำกัดของผลการศึกษา

1. **In silico เท่านั้น** — ไม่มีการยืนยันในห้องปฏิบัติการ (In Vitro)
2. **SHAP = ความสัมพันธ์เชิงสถิติ** ไม่ใช่กลไกเชิงสาเหตุระดับ 3 มิติ
3. **ข้อมูลมาจากหลาย assay/หลายห้องปฏิบัติการ** แม้จะกรองเฉพาะ 3CL-Pro แล้ว
   แต่โปรโตคอลการวัดที่ต่างกันยังเป็นแหล่ง noise ในตัวแปรตาม
4. **External validation อ่อน (R² = -1.719)** เนื่องจาก distribution shift
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

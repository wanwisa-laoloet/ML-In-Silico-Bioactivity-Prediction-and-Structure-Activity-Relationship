# แผนการดำเนินงานโครงงาน (Project Execution Plan)

**โครงงาน:** การประยุกต์ใช้การเรียนรู้ของเครื่องเพื่อทำนายฤทธิ์ทางชีวภาพและวิเคราะห์ความสัมพันธ์ระหว่างโครงสร้างกับฤทธิ์ของสารยับยั้ง SARS-CoV-2 Main Protease แบบ In Silico

**English:** Machine Learning Application for In Silico Bioactivity Prediction and Structure–Activity Relationship Analysis of SARS-CoV-2 Main Protease Inhibitors

| หัวข้อ | รายละเอียด |
|---|---|
| ผู้จัดทำ | 6633800639-0 นายรุ่งธรรม ล้อมในเมือง / 683380668-7 นางสาววรรณวิษา เลาเลิศ |
| อาจารย์ที่ปรึกษา | อ.ธนพล ตั้งชูพงศ์ |
| รายวิชา | CP352101 วิทยาการข้อมูลขั้นแนะนำ ภาคเรียนที่ 1 ปีการศึกษา 2569 |
| เอกสารอ้างอิงหลัก | `Docs/Machine Learning Application for In Silico Bioactivity Prediction and Structure-Activity Relationship.pdf` |
| วันที่จัดทำแผน | 8 กันยายน 2569 |

> ### สรุปย่อสำหรับผู้บริหาร (Executive Summary)
>
> ข้อเสนอโครงงานกำหนดให้ใช้ข้อมูล **SARS-CoV-2 Main Protease (Mpro/3CLpro) จาก ChEMBL ที่มี IC50 หน่วย nM** แต่จากการตรวจสอบชุดข้อมูลใน `Data_set/` พบว่า **ไม่มีไฟล์ใดตรงตามข้อกำหนดนี้**
>
> - `mergedreal.csv` / `xfinal.csv` มี 469 แถว แต่เป็นโมเลกุลซ้ำกันเพียง **10 ชนิด** (กรดอะมิโน/สารธรรมชาติ) → **ใช้ไม่ได้**
> - `DDH Data*.csv` เป็นสาร COVID-19 จริง 104 ตัว (ใช้ได้ 94) แต่หน่วยเป็น **µM ไม่ใช่ nM** และจำนวนน้อยเกินไปสำหรับ 1,028 features
> - `molecule_classification_dataset.csv` เป็นเป้าหมาย DRD2/HDAC1/BACE1 ไม่ใช่ COVID และไม่มีค่าฤทธิ์เชิงปริมาณ
>
> **ข้อเสนอแนะ:** ดึงข้อมูลหลักใหม่จาก ChEMBL target **`CHEMBL4523582`** (SARS-CoV-2 Replicase polyprotein 1ab / 3CL-Pro, UniProt P0DTD1) ซึ่งตรวจสอบผ่าน ChEMBL API แล้วเมื่อ 8 ก.ย. 2569 พบว่ามี **IC50 จำนวน 6,220 รายการ (หน่วย nM 6,050 รายการ, relation `=` 4,867 รายการ)** แล้วใช้ชุดข้อมูลเดิมเป็น External Validation และ Chemical-Space Reference ตามที่ระบุในส่วนที่ 4

---

## สารบัญ

1. [ส่วนที่ 1: สรุปสาระสำคัญจากข้อเสนอโครงงาน](#ส่วนที่-1-สรุปสาระสำคัญจากข้อเสนอโครงงาน)
2. [ส่วนที่ 2: ผลการตรวจสอบชุดข้อมูลที่มีอยู่ (Data Audit)](#ส่วนที่-2-ผลการตรวจสอบชุดข้อมูลที่มีอยู่-data-audit)
3. [ส่วนที่ 3: การวิเคราะห์ช่องว่าง (Gap Analysis)](#ส่วนที่-3-การวิเคราะห์ช่องว่าง-gap-analysis)
4. [ส่วนที่ 4: กลยุทธ์ข้อมูลที่เสนอ (Data Strategy)](#ส่วนที่-4-กลยุทธ์ข้อมูลที่เสนอ-data-strategy)
5. [ส่วนที่ 5: โครงสร้างโปรเจกต์และสภาพแวดล้อม](#ส่วนที่-5-โครงสร้างโปรเจกต์และสภาพแวดล้อม)
6. [ส่วนที่ 6: แผนการดำเนินงานรายเฟส (Phase 0–8)](#ส่วนที่-6-แผนการดำเนินงานรายเฟส-phase-08)
7. [ส่วนที่ 7: โปรโตคอลการสร้างแบบจำลองโดยละเอียด](#ส่วนที่-7-โปรโตคอลการสร้างแบบจำลองโดยละเอียด)
8. [ส่วนที่ 8: โปรโตคอล SHAP และการวิเคราะห์ SAR](#ส่วนที่-8-โปรโตคอล-shap-และการวิเคราะห์-sar)
9. [ส่วนที่ 9: การจับคู่วัตถุประสงค์กับผลผลิต (Traceability Matrix)](#ส่วนที่-9-การจับคู่วัตถุประสงค์กับผลผลิต-traceability-matrix)
10. [ส่วนที่ 10: ความเสี่ยงและแนวทางรับมือ](#ส่วนที่-10-ความเสี่ยงและแนวทางรับมือ)
11. [ส่วนที่ 11: ตารางเวลาดำเนินงาน](#ส่วนที่-11-ตารางเวลาดำเนินงาน)
12. [ภาคผนวก ก: คำสั่งและสคริปต์อ้างอิง](#ภาคผนวก-ก-คำสั่งและสคริปต์อ้างอิง)

---

## ส่วนที่ 1: สรุปสาระสำคัญจากข้อเสนอโครงงาน

### 1.1 วัตถุประสงค์ 3 ข้อ (ตามหัวข้อ 3 ของ PDF)

| # | วัตถุประสงค์ | ตัวชี้วัดความสำเร็จ (Definition of Done) |
|---|---|---|
| **O1** | พัฒนาและเปรียบเทียบประสิทธิภาพแบบจำลอง ML ในการทำนาย pIC50 ของสารยับยั้ง SARS-CoV-2 Mpro | มีตารางเปรียบเทียบ Dummy / RF / XGBoost / SVR บน Test Set ด้วย R², RMSE, MAE พร้อม Residual Plot |
| **O2** | สร้างกระบวนการเตรียมและสกัดคุณลักษณะแบบผสม (Hybrid Features) จาก Morgan FP + Lipinski Descriptors ผ่าน Scikit-learn Pipeline ที่ป้องกัน Data Leakage | มี `Pipeline(ColumnTransformer → Model)` ที่ scaling เกิดเฉพาะใน train fold, พิสูจน์ได้ด้วยการเทียบผลกับกรณี scale-ก่อน-split |
| **O3** | วิเคราะห์และอธิบายความสัมพันธ์ระหว่างคุณลักษณะโครงสร้างกับฤทธิ์ทางชีวภาพ (SAR) ด้วย SHAP + การวาดภาพโครงสร้างย่อย 2 มิติ | มี Bar plot Mean\|SHAP\| Top-20, Beeswarm plot, และภาพ `DrawMorganBit` ของบิตสำคัญ ≥ 6 บิต พร้อมคำอธิบายเชิงเคมี |

### 1.2 ข้อกำหนดทางเทคนิคที่ระบุไว้ในเอกสาร

| หัวข้อ | ข้อกำหนดจาก PDF |
|---|---|
| แหล่งข้อมูล | ChEMBL — SARS-CoV-2 Main Protease (Mpro / 3CLpro) |
| ตัวแปรตาม (Target) | `pIC50 = 9 − log10(IC50[nM]) = −log10(IC50[M])` |
| ประเภทปัญหา | Regression |
| Features | Morgan Fingerprints (ECFP4, r=2, **1,024 bits**) + Lipinski Descriptors 4 ตัว (MW, LogP, HBD, HBA) = **1,028 features** |
| Data Cleaning | กรอง IC50 หน่วย nM, Salt Stripping, Canonical SMILES, Deduplication ด้วยค่าเฉลี่ย |
| การแบ่งข้อมูล | Train 75% / Test 25% |
| Cross-Validation | 5-Fold CV บน Training Set |
| Models | Dummy Regressor (Baseline), Random Forest, XGBoost, SVR (RBF kernel) |
| Tuning | GridSearchCV บน Train Set (ระบุ Random Forest เป็นหลัก) |
| Metrics | R², RMSE, MAE + Residual Plot |
| XAI | SHAP → `Importance_j = (1/N) Σ|φ_j^(i)|` (Mean Absolute SHAP) |
| SAR | Bit Mapping ผ่าน `RDKit.DrawMorganBit` |
| Leakage Prevention | `ColumnTransformer` + `Pipeline` (StandardScaler เฉพาะ numeric descriptors) |
| Software | Python 3, RDKit, Scikit-learn, XGBoost, SHAP, Pandas, NumPy, Matplotlib, Seaborn |
| Hardware | CPU 4 Cores, RAM 16 GB |
| ผลผลิตที่คาดหวัง | `clean_bioactivity_master.csv`, ตารางผลโมเดล, ภาพ SAR |

### 1.3 ขอบเขตและข้อจำกัดที่ประกาศไว้ (หัวข้อ 6)

- **In Scope:** In silico เท่านั้น, 2D descriptors, Regression, 4 โมเดล, 5-Fold CV
- **Out of Scope:** In vitro validation, Molecular Docking, Molecular Dynamics, 3D descriptors
- **ข้อควรระวังที่เอกสารระบุเอง:** SHAP อธิบายได้เพียง **ความสัมพันธ์เชิงสถิติ (Statistical Association)** ไม่ใช่ **กลไกเชิงสาเหตุ (Causal Mechanism)** — ต้องคงข้อความนี้ไว้ในรายงานผล

---

## ส่วนที่ 2: ผลการตรวจสอบชุดข้อมูลที่มีอยู่ (Data Audit)

ตรวจสอบเมื่อ 8 กันยายน 2569 ด้วย pandas 2.3.2 / Python 3.13.5

### 2.1 ภาพรวมไฟล์ทั้งหมด

| # | ไฟล์ | ขนาด | แถว × คอลัมน์ | เป้าหมายชีวภาพ | มีค่าฤทธิ์? |
|---|---|---|---|---|---|
| D1 | `Bioactivity Dataset with Molecular Fingerprints/mergedreal.csv` | 3.9 MB | 469 × 4,264 | ไม่ระบุ (สารธรรมชาติ) | `pActivity` |
| D2 | `Bioactivity Dataset with Molecular Fingerprints/xfinal.csv` | 297 KB | 469 × 302 | ไม่ระบุ (เหมือน D1) | `pActivity` |
| D3 | `COVID-19 Drug Discovery Data/DDH Data.csv` | 7.3 KB | 104 × 3 | **COVID-19 (Mpro)** | `pIC50 (IC50 in microM)` |
| D4 | `COVID-19 Drug Discovery Data/DDH Data with Properties.csv` | 47 KB | 104 × 40 | **COVID-19 (Mpro)** | `pIC50` (เหมือน D3) |
| D5 | `MultiTarget Bioactivity ChEMBL/molecule_classification_dataset.csv` | 126 KB | 1,333 × 10 | DRD2 / HDAC1 / BACE1 | **ไม่มี** (มีแต่ label โปรตีน) |

### 2.2 D1 + D2 — Bioactivity Dataset with Molecular Fingerprints

**โครงสร้าง**

- `mergedreal.csv` = `SMILES`, `pActivity` + fingerprint 4,262 คอลัมน์:
  `FP_Morgan_0..2047` (2,048 bits) + `FP_MACCS_0..165` (166 bits) + `FP_RDKit_0..2047` (2,048 bits)
- `xfinal.csv` = `SMILES`, `pActivity` + fingerprint **300 คอลัมน์** (RDKit 257, Morgan 25, MACCS 18)
  - ตรวจสอบแล้ว: คอลัมน์ทั้ง 300 เป็น **subset ที่แท้จริงของ `mergedreal.csv`** และลำดับแถว SMILES ตรงกันทุกแถว
  - แปลว่า `xfinal.csv` คือผลลัพธ์ของ **Feature Selection ที่ทำบนข้อมูลทั้งชุดก่อนแบ่ง Train/Test**

**ปัญหาที่พบ (ร้ายแรง)**

| ปัญหา | รายละเอียด | ผลกระทบ |
|---|---|---|
| **โมเลกุลซ้ำมหาศาล** | 469 แถว แต่มี SMILES ไม่ซ้ำเพียง **10 ตัว** (แถวซ้ำ 459 แถว) | ขนาดตัวอย่างจริง N = 10 |
| **การกระจายไม่สมดุลรุนแรง** | 1 โมเลกุลปรากฏ 152 ครั้ง, อีกตัว 134 ครั้ง, อีกตัว 120 ครั้ง | ถ้าสุ่มแบ่ง train/test โมเลกุลเดียวกันจะอยู่ทั้งสองฝั่ง → **R² จะสูงปลอม ๆ เกือบ 1.00** |
| **ไม่ใช่สารยับยั้ง Mpro** | โมเลกุลคือ glycine, alanine, β-alanine, proline, serotonin, nicotine, linalool, eucalyptol, สาร indole/alkaloid | **ไม่ตรงวัตถุประสงค์ O1 เลย** |
| **Feature Selection Leakage ฝังมาแล้ว** | `xfinal.csv` เลือก feature จากข้อมูลเต็ม | **ขัดกับวัตถุประสงค์ O2 โดยตรง** |
| **Fingerprint spec ไม่ตรง** | Morgan 2,048 bits (PDF ระบุ 1,024 bits) และไม่มี Lipinski Descriptors | ต้องคำนวณใหม่อยู่ดี |

**รายละเอียดโมเลกุลทั้ง 10 ตัว**

| SMILES | จำนวนแถว | pActivity |
|---|---|---|
| `CC(=CCCC(C)(C=C)O)C` (Linalool) | 152 | 3.850 |
| `C1C[C@H](NC1)C(=O)O` (Proline) | 134 | 4.764 |
| `C(C(=O)O)N` (Glycine) | 120 | 5.905 |
| `C1=CC2=C(C=C1O)C(=CN2)CCN` (Serotonin) | 36 | 6.442 |
| `C1CN2CC3=CCO[C@H]4CC(=O)N5...` (Alkaloid) | 16 | 8.205 |
| `CN1CCC[C@H]1C2=CN=CC=C2` (Nicotine) | 6 | 4.149 |
| `C(CN)C(=O)O` (β-Alanine) | 2 | 4.936 |
| `CC1(C2CCC(O1)(CC2)C)C` (Eucalyptol) | 1 | 3.590 |
| `CC1=NC=CC2=C1NC3=C2C=CC(=C3)OC` (Harmine-like) | 1 | 5.000 |
| `C[C@@H](C(=O)O)N` (Alanine) | 1 | 3.703 |

> **คำตัดสิน (Verdict): ไม่ใช้เป็นชุดข้อมูลหลัก** ใช้ได้เพียงเป็น **Smoke-test fixture** สำหรับทดสอบว่าโค้ด Pipeline รันผ่าน (Phase 2) เท่านั้น และต้องบันทึกเหตุผลการตัดออกไว้ในรายงาน

### 2.3 D3 + D4 — COVID-19 Drug Discovery Data (DDH)

**โครงสร้าง**

- `DDH Data.csv`: `Compound No.`, `SMILES`, `pIC50 (IC50 in microM)` (+ คอลัมน์ว่าง `Unnamed: 3`)
- `DDH Data with Properties.csv`: SMILES เดียวกันทุกแถว + PubChem CID + descriptors 37 ตัว (MolecularWeight, XLogP, TPSA, HBondDonorCount, HBondAcceptorCount, RotatableBondCount, Complexity และ 3D features)

**สถิติสำคัญ**

| รายการ | ค่า |
|---|---|
| จำนวนสารทั้งหมด | 104 (SMILES ไม่ซ้ำครบ 104, ไม่มี salt/mixture `.`) |
| ค่า label ที่ใช้ได้ | **94** (อีก **10 แถวเป็นข้อความ `BLINDED`**) |
| ช่วง `pIC50 (µM basis)` | −2.699 ถึง 1.222 (mean −1.051, SD 0.888) |
| ช่วงหลังแปลงเป็น pIC50 มาตรฐาน (`+6`) | **3.301 ถึง 7.222** (mean 4.949) |
| เทียบเป็น IC50 | ~501 µM (อ่อนสุด) ถึง ~0.060 µM (แรงสุด) |
| ความยาว SMILES | 23–79 ตัวอักษร (mean 53.6) |
| ค่าหายใน D4 | 3 แถวไม่มี descriptor 2D เลย, 4 แถวไม่มี 3D features, **1 แถวมี `MolecularWeight = 0`** |

**ประเด็นสำคัญ — ความไม่ตรงกันของหน่วย**

```
PDF กำหนด:   pIC50 = 9 − log10(IC50[nM])   = −log10(IC50[M])
DDH ใช้:      pIC50 = −log10(IC50[µM])       = −log10(IC50[M]) − 6
∴ ต้องแปลง:   pIC50_standard = pIC50_DDH + 6
```

**ห้ามใช้ค่าดิบโดยไม่แปลง** มิฉะนั้นค่า target จะติดลบและตีความผิดทั้งฉบับ

> **คำตัดสิน (Verdict): ใช้ได้ แต่เป็นชุดรอง** เหมาะเป็น **External Validation Set** ที่ทรงคุณค่าที่สุด เพราะเป็นสาร COVID-19 จริงและมาจากแหล่งอิสระจาก ChEMBL (ต้องตรวจ overlap ด้วย InChIKey ก่อน) แต่ N = 94 น้อยเกินไปที่จะเทรนกับ 1,028 features ตามลำพัง

### 2.4 D5 — MultiTarget Bioactivity ChEMBL

| รายการ | ค่า |
|---|---|
| จำนวนแถว | 1,333 (ChEMBL ID ไม่ซ้ำ 1,333, SMILES ไม่ซ้ำ 1,331 → ซ้ำ 2) |
| คอลัมน์ | `ChEMBL ID`, `Name`, `SMILES`, `Molecular Weight`, `LogP`, `HBA`, `HBD`, `TPSA`, `Protein` |
| การกระจายเป้าหมาย | **DRD2 = 1,273 / HDAC1 = 58 / BACE1 = 2** (ไม่สมดุลรุนแรง) |
| ค่าหาย | `Name` หาย 1,206 แถว, descriptors หาย 8–10 แถว |
| ค่าฤทธิ์เชิงปริมาณ | **ไม่มี** — มีเพียง label ชื่อโปรตีน |
| SMILES ซ้อนทับกับ DDH | **0 แถว** |
| ช่วง MW | 135–1,659 (มี outlier > 1,000 Da ซึ่งเกินกฎ Lipinski มาก) |

> **คำตัดสิน (Verdict): ไม่ใช้กับ O1** เพราะเป็น classification คนละเป้าหมายและไม่มี IC50 แต่ **ใช้ได้ดีกับ 2 งานเสริม**: (ก) เป็นชุด Chemical Space เปรียบเทียบใน EDA ว่าสาร Mpro มีลักษณะต่างจากสารเป้าหมายอื่นอย่างไร (ข) เป็นชุดทดสอบ **Applicability Domain** ว่าโมเดลรู้ตัวหรือไม่เมื่อเจอสารนอกโดเมน

### 2.5 สรุปคะแนนความเหมาะสม

| ชุดข้อมูล | ตรงเป้าหมาย Mpro | มี IC50 nM | N เพียงพอ | ปลอด Leakage | **สรุป** |
|---|---|---|---|---|---|
| D1/D2 (Fingerprints) | ✗ | ✗ | ✗ (N=10) | ✗ | **ตัดออก** |
| D3/D4 (DDH COVID) | ✓ | ✗ (µM) | ✗ (N=94) | ✓ | **ใช้เป็น External Validation** |
| D5 (MultiTarget) | ✗ | ✗ | ✓ | ✓ | **ใช้เป็น Chemical Space / AD** |
| **ChEMBL4523582 (ต้องดึงใหม่)** | **✓** | **✓ (4,867 rec.)** | **✓** | **✓** | **ชุดข้อมูลหลัก** |

---

## ส่วนที่ 3: การวิเคราะห์ช่องว่าง (Gap Analysis)

| # | ข้อกำหนดใน PDF | สถานะปัจจุบัน | ระดับช่องว่าง | การแก้ไข |
|---|---|---|---|---|
| G1 | ข้อมูล SARS-CoV-2 Mpro จาก ChEMBL | ไม่มีไฟล์ใดตรง | **วิกฤต** | Phase 1: ดึงจาก ChEMBL API target `CHEMBL4523582` |
| G2 | IC50 หน่วย nM | DDH เป็น µM | ปานกลาง | แปลง `+6` และบันทึก provenance |
| G3 | Morgan FP 1,024 bits (r=2) | D1 มี 2,048 bits ที่ไม่ทราบ parameter | ปานกลาง | คำนวณใหม่ด้วย `GetMorganGenerator(radius=2, fpSize=1024)` |
| G4 | Lipinski Descriptors 4 ตัว | D4 มีจาก PubChem แต่มีค่าหาย/ค่า 0 | ปานกลาง | คำนวณใหม่ด้วย `RDKit.Descriptors` (MolWt, MolLogP, NumHDonors, NumHAcceptors) |
| G5 | ป้องกัน Data Leakage | `xfinal.csv` มี selection leakage ฝังมาแล้ว | **สูง** | ห้ามใช้ `xfinal.csv`; ทำ selection (ถ้ามี) ภายใน Pipeline เท่านั้น |
| G6 | Deduplication ด้วยค่าเฉลี่ย | D1 ซ้ำ 459 แถวโดยไม่ยุบ | **สูง** | ยุบด้วย median ต่อ InChIKey และตรวจ standard deviation |
| G7 | Library: RDKit, XGBoost, SHAP | **ยังไม่ติดตั้ง** (มีแค่ pandas 2.3.2, sklearn 1.7.2) | ปานกลาง | Phase 0: `pip install rdkit xgboost shap` |
| G8 | Train 75 / Test 25 | — | ต่ำ | ใช้ตามระบุ + เพิ่ม Scaffold Split เป็นการทดสอบเข้มงวด |
| G9 | ค่า `BLINDED` ใน DDH | 10 แถว | ต่ำ | คัดออกจาก supervised set; เก็บไว้เป็น prediction demo |

---

## ส่วนที่ 4: กลยุทธ์ข้อมูลที่เสนอ (Data Strategy)

ใช้สถาปัตยกรรม **3 ชั้น** เพื่อให้ทั้งตรงวัตถุประสงค์เดิมและใช้ประโยชน์จากไฟล์ที่มีอยู่ให้คุ้มค่า

```
┌─────────────────────────────────────────────────────────────────┐
│ ชั้นที่ 1 — PRIMARY (ใช้เทรน + ทดสอบ)                             │
│ ChEMBL CHEMBL4523582 : SARS-CoV-2 Replicase polyprotein 1ab      │
│ (3CL-Pro / Mpro), UniProt P0DTD1                                 │
│ IC50 ทั้งหมด 6,220 | nM 6,050 | relation "=" 4,867 รายการ         │
│ → หลัง cleaning + dedup คาดได้ ~1,500–2,500 โมเลกุลไม่ซ้ำ         │
│ → ส่งออกเป็น data/processed/clean_bioactivity_master.csv          │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ ชั้นที่ 2 — EXTERNAL VALIDATION (ไม่แตะจนกว่าจะ freeze โมเดล)      │
│ DDH Data.csv (94 สาร, แปลง +6) + descriptors จาก DDH Properties   │
│ → ทดสอบ generalization ข้ามแหล่งข้อมูล/ข้ามห้องปฏิบัติการ          │
│ → ตัดสาร InChIKey ที่ซ้ำกับชั้นที่ 1 ออกก่อนรายงานผล               │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ ชั้นที่ 3 — CONTEXT / APPLICABILITY DOMAIN                        │
│ molecule_classification_dataset.csv (DRD2/HDAC1/BACE1, 1,333)    │
│ → เปรียบเทียบ Chemical Space (MW vs LogP) ใน EDA                  │
│ → ทดสอบว่าโมเดลให้ความมั่นใจต่ำลงเมื่ออยู่นอก AD หรือไม่           │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│ ตัดออก — mergedreal.csv / xfinal.csv                             │
│ เหตุผล: N จริง = 10, ไม่ใช่ Mpro, มี selection leakage             │
│ ใช้ได้เพียง smoke-test ว่าโค้ด Pipeline รันไม่ error               │
└─────────────────────────────────────────────────────────────────┘
```

### 4.1 แผนสำรอง (Contingency Plan B)

หากอาจารย์ที่ปรึกษากำหนดว่า **ห้ามดึงข้อมูลใหม่** ให้ใช้ DDH 94 สารเป็นชุดหลัก โดยปรับเงื่อนไขดังนี้ และต้องระบุเป็น Limitation ในรายงาน:

- ลด Morgan FP จาก 1,024 → 256 bits และตัด bit ที่ variance = 0 ภายใน Pipeline (คาดเหลือ ~150–250 features)
- ใช้ **Repeated Nested CV (5×5, 20 รอบ)** แทนการแบ่ง 75/25 ครั้งเดียว เพราะ test set จะเหลือเพียง 24 สาร ซึ่งค่า R² จะผันผวนสูงมาก
- รายงานผลเป็น **ช่วงความเชื่อมั่น (mean ± SD ข้าม repeats)** ไม่ใช่ตัวเลขเดี่ยว
- เพิ่ม Ridge / ElasticNet เข้าไปในกลุ่มเปรียบเทียบ เพราะโมเดลเชิงเส้นที่มี regularization เหมาะกับ N น้อย

---

## ส่วนที่ 5: โครงสร้างโปรเจกต์และสภาพแวดล้อม

### 5.1 โครงสร้างโฟลเดอร์เป้าหมาย

```
SilicoBioactivityPredict/
├── Data_set/                              # ข้อมูลดิบเดิม (READ-ONLY ห้ามแก้)
│   ├── Bioactivity Dataset with Molecular Fingerprints/
│   ├── COVID-19 Drug Discovery Data/
│   └── MultiTarget Bioactivity ChEMBL/
├── data/
│   ├── raw/
│   │   └── chembl_CHEMBL4523582_ic50_raw.csv      # ผลดิบจาก API (Phase 1)
│   ├── interim/
│   │   ├── chembl_after_filter.csv
│   │   └── ddh_standardized.csv
│   └── processed/
│       ├── clean_bioactivity_master.csv           # ★ ผลผลิตตามข้อ 8.3 ของ PDF
│       ├── features_hybrid.parquet                # 1,028 features
│       ├── bitinfo.pkl                            # สำหรับ DrawMorganBit
│       ├── external_ddh_features.parquet
│       └── predictions_test.csv
├── notebooks/
│   ├── 00_data_audit.ipynb
│   ├── 01_data_acquisition.ipynb
│   ├── 02_cleaning_standardization.ipynb
│   ├── 03_feature_extraction.ipynb
│   ├── 04_eda.ipynb
│   ├── 05_modeling_cv.ipynb
│   ├── 06_tuning_final_eval.ipynb
│   ├── 07_shap_sar.ipynb
│   └── 08_external_validation.ipynb
├── src/
│   ├── config.py            # ค่าคงที่ทั้งหมด (RANDOM_STATE=42, FP_BITS=1024, ...)
│   ├── acquire.py           # ดึง ChEMBL
│   ├── clean.py             # standardize, salt strip, dedup, pIC50
│   ├── features.py          # Morgan FP + Lipinski
│   ├── pipelines.py         # ColumnTransformer + estimators
│   ├── evaluate.py          # metrics, residual plots
│   └── explain.py           # SHAP + DrawMorganBit
├── reports/
│   ├── figures/
│   │   ├── fig01_pic50_distribution.png
│   │   ├── fig02_chemical_space_mw_logp.png
│   │   ├── fig03_model_comparison.png
│   │   ├── fig04_residual_plot.png
│   │   ├── fig05_shap_bar_top20.png
│   │   ├── fig06_shap_beeswarm.png
│   │   └── fig07_bit_XXXX_substructure.png
│   ├── tables/
│   │   ├── tbl01_dataset_summary.csv
│   │   ├── tbl02_cv_results.csv
│   │   └── tbl03_test_results.csv
│   └── final_report.md
├── models/
│   └── best_model_<algo>_<date>.joblib
├── Docs/
│   ├── Machine Learning Application ... .pdf
│   └── Project_Plan.md                    # ← เอกสารนี้
├── requirements.txt
└── README.md
```

### 5.2 สภาพแวดล้อม — สถานะปัจจุบันและสิ่งที่ต้องติดตั้ง

**สถานะ: ติดตั้งเสร็จแล้วใน `.venv/` เฉพาะโปรเจกต์ (Python 3.13.5)** — ตรึงเวอร์ชันไว้ใน `requirements.txt`

| Library | ต้องการ | เวอร์ชันที่ติดตั้งจริงใน `.venv` | ใช้ทำอะไร |
|---|---|---|---|
| Python | 3.x | **3.13.5** | — |
| pandas | ✓ | **2.3.2** | จัดการข้อมูลตาราง |
| numpy | ✓ | **2.2.6** | คำนวณเชิงตัวเลข |
| scikit-learn | ✓ | **1.7.2** | Pipeline, CV, RF, SVR, metrics |
| **RDKit** | ✓ | **2026.3.6** | Morgan FP, Lipinski, DrawMorganBit |
| **XGBoost** | ✓ | **3.4.1** | XGBoost Regressor |
| **SHAP** | ✓ | **0.52.0** | TreeExplainer สำหรับ SAR |
| matplotlib | ✓ | **3.11.1** | สร้างกราฟ |
| seaborn | ✓ | **0.13.2** | ธีมและ statistical plots |
| scipy | (เพิ่ม) | **1.16.3** | Mann-Whitney U, Spearman ρ |
| pyarrow | (เพิ่ม) | **25.0.1** | อ่าน/เขียน parquet |
| joblib | (เพิ่ม) | **1.5.2** | บันทึก/โหลดโมเดล |

> **ผลจริง (ความเสี่ยง R6 ไม่เกิดขึ้น):** `rdkit`, `xgboost` และ `shap` ติดตั้งบน **Python 3.13.5 ได้สำเร็จทั้งหมด** จึงไม่จำเป็นต้องถอยไปใช้ Python 3.11 ตามแผนสำรอง
>
> **หมายเหตุ:** ไม่ได้ติดตั้ง `chembl-webresource-client` เพราะ `src/acquire.py` เรียก ChEMBL REST API ผ่าน `urllib` ในไลบรารีมาตรฐานโดยตรง ซึ่งควบคุม paging และ retry ได้ละเอียดกว่า

**คำสั่งสร้าง venv ใหม่ (กรณีย้ายเครื่อง):**

```bash
python -m venv .venv
```

```bash
./.venv/Scripts/python.exe -m pip install -r requirements.txt
```

### 5.3 การควบคุมความสามารถทำซ้ำได้ (Reproducibility)

- กำหนด `RANDOM_STATE = 42` ที่ `src/config.py` และส่งต่อไปยังทุกจุดที่มีการสุ่ม (`train_test_split`, `KFold`, `RandomForestRegressor`, `XGBRegressor`)
- ตรึงเวอร์ชันใน `requirements.txt` ด้วย `pip freeze`
- บันทึก ChEMBL release และวันที่ดึงข้อมูลไว้ใน metadata ของ `clean_bioactivity_master.csv`
- เริ่มใช้ git (`git init`) — ปัจจุบันโฟลเดอร์ยังไม่เป็น repository

---

## ส่วนที่ 6: แผนการดำเนินงานรายเฟส (Phase 0–8)

### Phase 0 — เตรียมสภาพแวดล้อมและตรวจรับข้อมูล (2–3 วัน)

| งาน | รายละเอียด | ผลผลิต |
|---|---|---|
| T0.1 | สร้าง venv (แนะนำ Python 3.11) และติดตั้ง rdkit, xgboost, shap, matplotlib, seaborn | `requirements.txt` |
| T0.2 | `git init` + `.gitignore` (ยกเว้น `data/raw/`, `models/`) | repo พร้อมใช้ |
| T0.3 | สร้างโครงสร้างโฟลเดอร์ตาม 5.1 | โครงสร้างครบ |
| T0.4 | เขียน `notebooks/00_data_audit.ipynb` บันทึกผลตรวจสอบในส่วนที่ 2 อย่างเป็นทางการ | `tbl01_dataset_summary.csv` |
| T0.5 | **นำเสนออาจารย์ที่ปรึกษา**: ชี้แจงว่าชุดข้อมูลเดิมไม่ตรงวัตถุประสงค์ ขออนุมัติดึง ChEMBL ใหม่ | บันทึกการอนุมัติ |

**เกณฑ์ผ่าน:** `import rdkit, xgboost, shap` สำเร็จ + ได้รับอนุมัติแผนข้อมูล

---

### Phase 1 — การรวบรวมข้อมูล (Data Acquisition) (3–4 วัน)

**สอดคล้องกับ:** ขั้นตอนที่ 1 ของหัวข้อ 5.2 ใน PDF

| งาน | รายละเอียด |
|---|---|
| T1.1 | ยืนยัน target: `CHEMBL4523582` = *Replicase polyprotein 1ab*, SARS-CoV-2, UniProt **P0DTD1**, target_type = SINGLE PROTEIN |
| T1.2 | ดึง activities ด้วยเงื่อนไข: `target_chembl_id=CHEMBL4523582`, `standard_type=IC50`, `standard_units=nM`, `standard_relation="="` (คาด **4,867 รายการ**) |
| T1.3 | ดึงฟิลด์: `molecule_chembl_id`, `canonical_smiles`, `standard_value`, `standard_units`, `standard_relation`, `assay_chembl_id`, `assay_description`, `document_chembl_id`, `pchembl_value`, `data_validity_comment` |
| T1.4 | **กรอง assay เฉพาะ 3CL-Pro/Mpro**: ใช้ regex บน `assay_description` จับคำ `3CL`, `3C-like`, `Mpro`, `main protease` เพื่อตัด assay ของโดเมนอื่นใน polyprotein 1ab (เช่น PLpro, RdRp, helicase) ออก — **ขั้นตอนนี้สำคัญมาก** เพราะ CHEMBL4523582 ครอบคลุมทั้ง polyprotein |
| T1.5 | บันทึกดิบเป็น `data/raw/chembl_CHEMBL4523582_ic50_raw.csv` พร้อม timestamp |
| T1.6 | เตรียม DDH: อ่าน `DDH Data.csv`, ตัด `Unnamed: 3`, แยก 10 แถว `BLINDED` ออกเป็น `ddh_blinded.csv`, แปลง `pIC50 + 6` |

**เกณฑ์ผ่าน:** ได้ CSV ดิบ ≥ 3,000 แถวหลังกรอง assay + ไฟล์ DDH 94 แถวที่แปลงหน่วยแล้ว

---

### Phase 2 — การทำความสะอาดข้อมูล (Data Cleaning) (4–5 วัน)

**สอดคล้องกับ:** ขั้นตอนที่ 2 ของหัวข้อ 5.2 ใน PDF

ลำดับขั้นตอนที่ต้องทำ **ตามลำดับนี้เท่านั้น** (ลำดับสลับจะทำให้ dedup ผิด):

```
 1. ตัดแถวที่ canonical_smiles หรือ standard_value เป็น null
 2. ตัดแถวที่ data_validity_comment ไม่ว่าง (เช่น "Outside typical range", "Potential author error")
 3. ตัด standard_value <= 0
 4. RDKit: MolFromSmiles → ตัดแถวที่ parse ไม่ได้ (บันทึกจำนวนไว้)
 5. Salt Stripping: rdMolStandardize.LargestFragmentChooser
 6. Normalize + Uncharge: rdMolStandardize.Normalizer() แล้ว Uncharger()
 7. Canonical SMILES: Chem.MolToSmiles(mol, canonical=True)
 8. คำนวณ InChIKey เป็นกุญแจสำหรับ dedup (แม่นกว่า SMILES)
 9. คำนวณ pIC50 = 9 − log10(standard_value[nM])
10. Deduplication: group by InChIKey
      - ถ้า SD ของ pIC50 ในกลุ่ม > 1.0 log unit → ตัดทั้งกลุ่มทิ้ง (ข้อมูลขัดแย้งกันเอง)
      - ถ้า SD <= 1.0 → ใช้ค่า median (ทนต่อ outlier กว่า mean)
      - บันทึก n_records และ pIC50_std เป็นคอลัมน์เพื่อความโปร่งใส
11. กรองช่วงที่สมเหตุสมผล: 3.0 <= pIC50 <= 11.0
12. ตัดโมเลกุลที่มีอะตอมหนัก < 5 หรือ > 70
13. เพิ่ม bioactivity_class: active (pIC50 >= 6), intermediate (5–6), inactive (< 5) — ใช้เพื่อ EDA และ stratify
```

**ผลผลิต:** `data/processed/clean_bioactivity_master.csv` — ตรงตามข้อ 8.3 ของ PDF

| คอลัมน์ | ชนิด | คำอธิบาย |
|---|---|---|
| `molecule_chembl_id` | str | รหัส ChEMBL |
| `canonical_smiles` | str | SMILES มาตรฐานหลังทำความสะอาด |
| `inchikey` | str | กุญแจ dedup |
| `ic50_nM` | float | ค่ามัธยฐานหลังยุบซ้ำ |
| `pIC50` | float | `9 − log10(ic50_nM)` |
| `n_records` | int | จำนวน record ต้นทางที่ยุบมา |
| `pIC50_std` | float | ส่วนเบี่ยงเบนมาตรฐานภายในกลุ่ม |
| `bioactivity_class` | str | active / intermediate / inactive |
| `source` | str | `ChEMBL` หรือ `DDH` |

**T2.x เพิ่มเติม:** รัน pipeline เดียวกันกับ DDH 94 สาร แล้วตรวจ **overlap ด้วย InChIKey** — สารที่ซ้ำกับชุด ChEMBL ต้องตัดออกจาก external set ก่อนรายงานผล

**เกณฑ์ผ่าน:** ไม่มี InChIKey ซ้ำใน master, ไม่มี null ใน `pIC50`, มี log บันทึกจำนวนแถวที่หายไปในทุกขั้นตอน

---

### Phase 3 — การสกัดคุณลักษณะ (Feature Extraction) (2–3 วัน)

**สอดคล้องกับ:** ขั้นตอนที่ 3 ของหัวข้อ 5.2 และวัตถุประสงค์ **O2**

| งาน | รายละเอียด |
|---|---|
| T3.1 | **Morgan Fingerprints**: `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)` → `GetFingerprintAsNumPy()` ได้ 1,024 คอลัมน์ `FP_0..FP_1023` (ค่า 0/1) |
| T3.2 | เก็บ `bitInfoMap` ของทุกโมเลกุลไว้ใน `data/processed/bitinfo.pkl` — **จำเป็นสำหรับ Phase 7** (`DrawMorganBit` ต้องใช้) |
| T3.3 | **Lipinski Descriptors 4 ตัว** ด้วย RDKit: `Descriptors.MolWt`, `Descriptors.MolLogP` (Crippen), `Lipinski.NumHDonors`, `Lipinski.NumHAcceptors` → `MW`, `LogP`, `HBD`, `HBA` |
| T3.4 | ตรวจสอบว่า `MW` ไม่มีค่า 0 (ปัญหาที่พบใน DDH Properties) |
| T3.5 | รวมเป็นเมทริกซ์ **1,028 features** บันทึกเป็น parquet (เร็วและเล็กกว่า CSV มาก) |
| T3.6 | คำนวณ `lipinski_violations` (0–4) ไว้ใช้ใน EDA เท่านั้น — **ห้ามใส่เป็น feature** เพราะเป็นฟังก์ชันของ 4 ตัวแรก (redundant) |

> **ข้อควรระวังสำคัญ:** **ห้าม** ตัด fingerprint bit ที่ variance ต่ำนอก Pipeline บนข้อมูลทั้งชุด — นั่นคือ leakage แบบเดียวกับที่พบใน `xfinal.csv` หากต้องการทำ ให้ใส่ `VarianceThreshold` เป็น **step แรกใน Pipeline** เพื่อให้ fit เฉพาะ train fold

**เกณฑ์ผ่าน:** shape = (N, 1028), ไม่มี NaN, fingerprint เป็น 0/1 ล้วน, `bitinfo.pkl` เปิดใช้ได้

---

### Phase 4 — การสำรวจข้อมูล (EDA) (3 วัน)

**สอดคล้องกับ:** ขั้นตอนที่ 4 ของหัวข้อ 5.2

| งาน | รายละเอียด | ผลผลิต |
|---|---|---|
| T4.1 | Histogram + KDE ของ pIC50 พร้อมสถิติ (mean, median, SD, skewness) | `fig01_pic50_distribution.png` |
| T4.2 | Chemical Space Scatter: MW (แกน x) vs LogP (แกน y) ระบายสีตาม `bioactivity_class` + เส้นกรอบ Lipinski (MW=500, LogP=5) | `fig02_chemical_space_mw_logp.png` |
| T4.3 | Box plot เปรียบเทียบ MW, LogP, HBD, HBA ระหว่าง active / intermediate / inactive + Mann-Whitney U test | `fig02b_lipinski_by_class.png` |
| T4.4 | สถิติ Fingerprint: จำนวนบิตที่ variance = 0, การกระจายของ bit density ต่อโมเลกุล | ตารางในโน้ตบุ๊ก |
| T4.5 | **เปรียบเทียบ Chemical Space 3 ชุด**: ChEMBL Mpro vs DDH vs MultiTarget (DRD2/HDAC1) บนแกน MW–LogP เดียวกัน | `fig02c_chemspace_3sets.png` |
| T4.6 | หา Murcko scaffold ที่พบบ่อย 10 อันดับแรก และตรวจ outlier ด้วย Tanimoto similarity | ตาราง |

**เกณฑ์ผ่าน:** ครบ 5 ภาพ พร้อมคำบรรยายเชิงเคมี ≥ 3 บรรทัดต่อภาพ

---

### Phase 5 — การสร้าง Pipeline และ Cross-Validation (4–5 วัน)

**สอดคล้องกับ:** ขั้นตอนที่ 5 ของหัวข้อ 5.2 และวัตถุประสงค์ **O1, O2**

#### T5.1 การแบ่งข้อมูล

```python
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42,
    stratify=df['bioactivity_class']      # ทำให้สัดส่วน active/inactive เท่ากันสองฝั่ง
)
```

**เพิ่มเติม (แนะนำอย่างยิ่ง):** ทำ **Scaffold Split** (Bemis–Murcko) เป็นการทดสอบชุดที่สอง เพราะ random split ทำให้สารในอนุกรมเคมีเดียวกัน (analog series) กระจายทั้งสองฝั่ง → R² สูงเกินจริง Scaffold split สะท้อนสถานการณ์การทำนายสารใหม่ได้ตรงกว่า **รายงานทั้งสองค่าคู่กัน**

#### T5.2 โครงสร้าง Pipeline (หัวใจของวัตถุประสงค์ O2)

```python
preprocessor = ColumnTransformer(
    transformers=[
        ('fp',   'passthrough',    fp_cols),        # 1,024 bits — ไม่ต้อง scale (0/1 อยู่แล้ว)
        ('desc', StandardScaler(), ['MW', 'LogP', 'HBD', 'HBA']),
    ],
    remainder='drop'
)

pipe = Pipeline([
    ('prep',  preprocessor),
    ('model', estimator),
])
```

**เหตุผลที่ป้องกัน Leakage ได้:** `StandardScaler` ถูก `fit` ใหม่ทุกครั้งที่ `cross_val_score` เรียก `pipe.fit(X_train_fold)` ค่า mean/SD จึงคำนวณจาก **train fold เท่านั้น** ไม่รั่วไปยัง validation fold

#### T5.3 การทดลองพิสูจน์ผลของ Leakage (สำหรับ O2 โดยเฉพาะ)

รัน 2 กรณีเทียบกัน แล้วรายงานส่วนต่าง — เป็นหลักฐานเชิงประจักษ์ว่า Pipeline มีความจำเป็นจริง:

| กรณี | วิธีทำ | คาดการณ์ |
|---|---|---|
| A (ถูกต้อง) | Scale ภายใน Pipeline | R² ต่ำกว่าเล็กน้อยแต่ตรงความจริง |
| B (ผิด) | `StandardScaler().fit_transform(X)` ทั้งชุดก่อน split | R² สูงกว่าอย่างมีอคติ |

#### T5.4 โมเดล 4 ตัวและค่าเริ่มต้น

| โมเดล | Class | Baseline params |
|---|---|---|
| Dummy | `DummyRegressor(strategy='mean')` | — |
| Random Forest | `RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)` | — |
| XGBoost | `XGBRegressor(n_estimators=500, learning_rate=0.05, max_depth=6, subsample=0.8, colsample_bytree=0.8, random_state=42, n_jobs=4)` | — |
| SVR | `SVR(kernel='rbf', C=10, gamma='scale', epsilon=0.1)` | — |

#### T5.5 Cross-Validation

```python
cv = KFold(n_splits=5, shuffle=True, random_state=42)
scoring = {'r2': 'r2',
           'rmse': 'neg_root_mean_squared_error',
           'mae': 'neg_mean_absolute_error'}
cross_validate(pipe, X_train, y_train, cv=cv, scoring=scoring, return_train_score=True)
```

รายงาน **mean ± SD** ของทุก metric ทั้ง train และ validation (ส่วนต่าง train−val คือดัชนี overfitting)

**เกณฑ์ผ่าน:** ทั้ง 4 โมเดลรันครบ, `tbl02_cv_results.csv` สมบูรณ์, ทุกโมเดลที่ไม่ใช่ Dummy ต้องมี R² สูงกว่า Dummy อย่างชัดเจน (Dummy จะได้ R² ≈ 0 หรือติดลบเล็กน้อย)

---

### Phase 6 — Hyperparameter Tuning และการประเมินขั้นสุดท้าย (4 วัน)

**สอดคล้องกับ:** ขั้นตอนที่ 5–6 ของหัวข้อ 5.2

#### T6.1 GridSearchCV (บน Training Set เท่านั้น)

**Random Forest** (โมเดลหลักตามที่ PDF ระบุ):

```python
param_grid_rf = {
    'model__n_estimators':      [300, 500, 800],
    'model__max_depth':         [None, 15, 25],
    'model__min_samples_split': [2, 5, 10],
    'model__min_samples_leaf':  [1, 2, 4],
    'model__max_features':      ['sqrt', 0.3],
}
# 3×3×3×3×2 = 162 combos × 5 folds = 810 fits
```

**XGBoost:**

```python
param_grid_xgb = {
    'model__n_estimators':     [300, 600],
    'model__learning_rate':    [0.03, 0.05, 0.1],
    'model__max_depth':        [4, 6, 8],
    'model__subsample':        [0.8, 1.0],
    'model__colsample_bytree': [0.6, 0.8],
}
```

**SVR:**

```python
param_grid_svr = {
    'model__C':       [1, 10, 100],
    'model__gamma':   ['scale', 0.001, 0.01],
    'model__epsilon': [0.05, 0.1, 0.2],
}
```

> **ข้อจำกัดฮาร์ดแวร์ (4 cores / 16 GB):** ถ้า GridSearchCV ของ RF ใช้เวลาเกิน ~2 ชั่วโมง ให้เปลี่ยนเป็น `RandomizedSearchCV(n_iter=60)` ผลลัพธ์ใกล้เคียงกันมากในทางปฏิบัติ และบันทึกเหตุผลไว้ในรายงาน

#### T6.2 การประเมินบน Test Set (**แตะเพียงครั้งเดียว**)

| Metric | สูตร | ความหมาย |
|---|---|---|
| R² | `1 − SS_res/SS_tot` | สัดส่วนความแปรปรวนที่อธิบายได้ |
| RMSE | `sqrt(mean((y−ŷ)²))` | ความคลาดเคลื่อนหน่วย log |
| MAE | `mean(\|y−ŷ\|)` | ความคลาดเคลื่อนเฉลี่ยสัมบูรณ์ |
| **เพิ่ม** Spearman ρ | — | ความสามารถในการ **จัดอันดับ** สาร ซึ่งสำคัญกว่าค่าสัมบูรณ์ในงานคัดกรองยา |

#### T6.3 กราฟผลลัพธ์

- `fig03_model_comparison.png` — Bar chart R²/RMSE/MAE ของ 4 โมเดล (CV vs Test เคียงกัน)
- `fig04_residual_plot.png` — 2 แผง: (ก) Predicted vs Actual + เส้น y=x (ข) Residual vs Predicted + เส้น 0 พร้อมแถบ ±1 log unit
- `fig04b_residual_hist.png` — Histogram ของ residual ตรวจ normality

#### T6.4 การตีความเกณฑ์

| ช่วง Test R² | การตีความในบริบท QSAR |
|---|---|
| < 0.40 | อ่อน — ต้องอภิปรายสาเหตุ (ข้อมูลผสมหลาย assay, activity cliff) |
| 0.40–0.60 | ยอมรับได้สำหรับ QSAR ข้ามชุดข้อมูล ChEMBL |
| 0.60–0.75 | **ดี** — ระดับที่คาดหวังจากโครงงานนี้ |
| > 0.75 | ดีมาก — **ต้องตรวจ leakage ซ้ำ** ว่าไม่ได้มาจาก analog series ที่กระจายสองฝั่ง |

**เกณฑ์ผ่าน:** บันทึกโมเดลที่ดีที่สุดเป็น `.joblib`, มี `tbl03_test_results.csv`, มี `predictions_test.csv`

---

### Phase 7 — SHAP และการวิเคราะห์ SAR (5–6 วัน)

**สอดคล้องกับ:** ขั้นตอนที่ 7 ของหัวข้อ 5.2 และวัตถุประสงค์ **O3** — รายละเอียดเต็มอยู่ในส่วนที่ 8

---

### Phase 8 — External Validation, รายงาน และการส่งมอบ (5 วัน)

| งาน | รายละเอียด |
|---|---|
| T8.1 | ทำนาย DDH 94 สาร (ที่ไม่ซ้ำกับ training) ด้วยโมเดลที่ freeze แล้ว → รายงาน R², RMSE, Spearman ρ |
| T8.2 | อภิปราย: ถ้าผลตกลงมาก ให้วิเคราะห์ว่าเป็นเพราะ chemical space ต่างกัน หรือ assay protocol ต่างกัน (DDH อยู่ในช่วง µM ซึ่งอ่อนกว่าค่าเฉลี่ยของ ChEMBL) |
| T8.3 | ทำนาย 10 สาร `BLINDED` เป็น demo การใช้งานจริง (ไม่มี ground truth จึงรายงานเป็นค่าทำนายเท่านั้น) |
| T8.4 | ทดสอบ Applicability Domain กับ MultiTarget 1,333 สาร — แสดงว่า RF prediction variance สูงขึ้นเมื่ออยู่นอกโดเมน |
| T8.5 | เขียน `reports/final_report.md` ตามโครง: บทนำ → วิธีการ → ผล → อภิปราย → ข้อจำกัด → สรุป |
| T8.6 | จัดทำสไลด์นำเสนอ + ตรวจสอบว่าผลผลิตครบตามข้อ 8 ของ PDF ทั้ง 3 ข้อ |
| T8.7 | ตรวจ Reproducibility: ล้าง environment แล้วรันซ้ำจากศูนย์ ตัวเลขต้องตรงกันทุกตำแหน่ง |

---

## ส่วนที่ 7: โปรโตคอลการสร้างแบบจำลองโดยละเอียด

### 7.1 ตารางสรุปการทดลองทั้งหมด

| Exp. ID | ชุดข้อมูล | Split | Features | โมเดล | จุดประสงค์ |
|---|---|---|---|---|---|
| E1 | ChEMBL Mpro | Random 75/25 | 1,028 hybrid | 4 โมเดล | **ผลหลักตอบ O1** |
| E2 | ChEMBL Mpro | Scaffold 75/25 | 1,028 hybrid | 4 โมเดล | ทดสอบเข้มงวด (generalization จริง) |
| E3 | ChEMBL Mpro | Random 75/25 | FP 1,024 อย่างเดียว | RF | Ablation: Lipinski ช่วยหรือไม่ |
| E4 | ChEMBL Mpro | Random 75/25 | Lipinski 4 อย่างเดียว | RF | Ablation: FP ช่วยแค่ไหน |
| E5 | ChEMBL Mpro | Random 75/25 | 1,028 hybrid | RF (scale ก่อน split) | **พิสูจน์ leakage ตอบ O2** |
| E6 | DDH (external) | — | 1,028 hybrid | โมเดลที่ดีที่สุดจาก E1 | ทดสอบข้ามชุดข้อมูล |

> E3 และ E4 คือ **Ablation Study** ที่ตอบคำถามสำคัญว่า "Hybrid Features คุ้มค่าจริงหรือไม่" ซึ่งเป็นเหตุผลตั้งต้นของวัตถุประสงค์ O2

### 7.2 รูปแบบตารางผลลัพธ์ที่ต้องส่งมอบ

**ตารางที่ 2 — ผล 5-Fold CV บน Training Set**

| Model | CV R² (mean±SD) | CV RMSE | CV MAE | Train R² | Δ(Train−CV) |
|---|---|---|---|---|---|
| Dummy Regressor | | | | | |
| Random Forest | | | | | |
| XGBoost | | | | | |
| SVR (RBF) | | | | | |

**ตารางที่ 3 — ผลบน Test Set (25%)**

| Model | Test R² | RMSE | MAE | Spearman ρ | เวลาเทรน (s) |
|---|---|---|---|---|---|
| Dummy Regressor | | | | | |
| Random Forest (tuned) | | | | | |
| XGBoost (tuned) | | | | | |
| SVR (tuned) | | | | | |

---

## ส่วนที่ 8: โปรโตคอล SHAP และการวิเคราะห์ SAR

**ตอบวัตถุประสงค์ O3 โดยตรง**

### 8.1 การคำนวณ SHAP

| ขั้น | รายละเอียด |
|---|---|
| S1 | ใช้ `shap.TreeExplainer` กับโมเดล tree ที่ดีที่สุด (RF หรือ XGBoost) — เป็น exact algorithm สำหรับ tree ไม่ต้องประมาณ |
| S2 | คำนวณบน **Test Set** ตามที่ PDF ระบุ (ไม่ใช่ train — เพื่อวัดว่าโมเดลใช้ feature อะไรตอนเจอข้อมูลใหม่) |
| S3 | ถ้า Test Set ใหญ่จนช้า ให้สุ่มตัวอย่าง 300–500 แถวและระบุไว้ในรายงาน |
| S4 | คำนวณ `Importance_j = (1/N) Σ_i \|φ_j^(i)\|` ตามสูตรในหัวข้อ 4.5 ของ PDF |
| S5 | **หมายเหตุสำคัญ:** ต้องดึง feature matrix **หลังผ่าน `preprocessor`** ออกมาก่อน แล้วส่งเข้า TreeExplainer พร้อม feature names ที่ตรงลำดับ มิฉะนั้นชื่อบิตจะเพี้ยน |

### 8.2 ภาพผลลัพธ์ที่ต้องมี

| ภาพ | เนื้อหา |
|---|---|
| `fig05_shap_bar_top20.png` | Bar chart Mean\|SHAP\| 20 features สูงสุด แยกสี fingerprint bits vs Lipinski descriptors |
| `fig06_shap_beeswarm.png` | Beeswarm plot แสดงทิศทาง (บิต = 1 ทำให้ pIC50 เพิ่มหรือลด) |
| `fig06b_shap_dependence_logp.png` | Dependence plot ของ LogP — ตรวจว่ามีช่วง optimum หรือไม่ |

### 8.3 Bit Mapping → ภาพโครงสร้างย่อย 2 มิติ

ขั้นตอนสำหรับแต่ละบิตสำคัญ (ทำอย่างน้อย **6 บิต**: 3 บิตที่ผลักดัน pIC50 ขึ้น, 3 บิตที่ผลักลง):

```
1. เลือก bit index j จาก Top-20 Mean|SHAP| ที่เป็น fingerprint
2. ค้นหาโมเลกุลใน test set ที่มีบิตนั้น = 1 และมี |SHAP| สูงสุด
3. ดึง bitInfo ของโมเลกุลนั้นจาก data/processed/bitinfo.pkl
4. วาดด้วย Draw.DrawMorganBit(mol, j, bitInfo)
5. บันทึกเป็น reports/figures/fig07_bit_<j>_substructure.png
6. เขียนคำอธิบายเชิงเคมี: หมู่ฟังก์ชันคืออะไร, พบใน scaffold ใดบ้าง,
   กี่ % ของสาร active มีบิตนี้เทียบกับสาร inactive
```

### 8.4 ตารางสรุป SAR ที่ต้องส่งมอบ

| Bit | Mean\|SHAP\| | ทิศทาง | โครงสร้างย่อย | % ในกลุ่ม Active | % ในกลุ่ม Inactive | การตีความเชิงเคมี |
|---|---|---|---|---|---|---|
| FP_xxx | | ↑ เพิ่มฤทธิ์ | (ภาพ) | | | |
| FP_yyy | | ↓ ลดฤทธิ์ | (ภาพ) | | | |

### 8.5 ข้อความบังคับในรายงาน (ตามข้อจำกัด 6.2.2 ของ PDF)

> "ผลการวิเคราะห์ SHAP ในงานนี้แสดง **ความสัมพันธ์เชิงสถิติ (Statistical Association)** ระหว่างชิ้นส่วนโมเลกุลกับค่าที่แบบจำลองทำนาย ไม่สามารถอ้างเป็น **กลไกการจับทางชีวภาพเชิงสาเหตุ (Causal Binding Mechanism)** ในระดับ 3 มิติได้ การยืนยันกลไกที่แท้จริงจำเป็นต้องอาศัย Molecular Docking, Molecular Dynamics และการทดสอบ In Vitro ซึ่งอยู่นอกขอบเขตของโครงงานนี้"

---

## ส่วนที่ 9: การจับคู่วัตถุประสงค์กับผลผลิต (Traceability Matrix)

| วัตถุประสงค์ | Phase ที่เกี่ยวข้อง | ผลผลิตที่พิสูจน์ความสำเร็จ |
|---|---|---|
| **O1** — พัฒนาและเปรียบเทียบโมเดล ML ทำนาย pIC50 | 1, 2, 5, 6 | `tbl02_cv_results.csv`, `tbl03_test_results.csv`, `fig03_model_comparison.png`, `fig04_residual_plot.png`, `models/best_model_*.joblib` |
| **O2** — สร้างกระบวนการ Hybrid Features ผ่าน Pipeline ที่ปลอด Leakage | 3, 5 (T5.2, T5.3), 7 (E3–E5) | `src/pipelines.py`, ผลเปรียบเทียบ E1 vs E5 (leakage demo), ผล Ablation E3/E4, `features_hybrid.parquet` |
| **O3** — วิเคราะห์ SAR ด้วย SHAP + ภาพโครงสร้าง 2 มิติ | 7 (Phase 7 ทั้งหมด) | `fig05_shap_bar_top20.png`, `fig06_shap_beeswarm.png`, `fig07_bit_*.png` ≥ 6 ภาพ, ตาราง SAR ในหัวข้อ 8.4 |
| **ประโยชน์ที่คาดว่าจะได้รับ ข้อ 3** — ชุดข้อมูลมาตรฐาน | 2 | `data/processed/clean_bioactivity_master.csv` |

---

## ส่วนที่ 10: ความเสี่ยงและแนวทางรับมือ

| ID | ความเสี่ยง | ระดับ | ผลกระทบ | แนวทางรับมือ |
|---|---|---|---|---|
| **R1** | ชุดข้อมูลที่มีอยู่ไม่ตรงวัตถุประสงค์เลย | **สูงมาก** | ทำ O1 ไม่ได้ | ดึง ChEMBL `CHEMBL4523582` ใหม่ (Phase 1) — ยืนยันแล้วว่ามี 4,867 records พร้อมใช้ |
| **R2** | assay ใน CHEMBL4523582 ปนหลายโดเมนของ polyprotein 1ab (PLpro, RdRp) | สูง | โมเดลเรียนรู้จากเป้าหมายผสม | กรอง `assay_description` ด้วย regex `3CL|3C-like|Mpro|main protease` (T1.4) และรายงานจำนวนที่คัดออก |
| **R3** | Activity Cliff — โมเลกุลคล้ายกันมากแต่ฤทธิ์ต่างกันหลาย log | สูง | RMSE สูงโดยไม่ใช่ความผิดของโมเดล | ระบุคู่ activity cliff ด้วย Tanimoto > 0.85 และ Δ pIC50 > 1.5 แล้วอภิปรายแยกใน Residual analysis |
| **R4** | Random split ทำให้ analog series กระจายสองฝั่ง → R² สูงเกินจริง | สูง | ผลลวง | ทำ Scaffold Split คู่ขนาน (E2) และรายงานทั้งสองค่า |
| **R5** | ข้อมูลจากหลายห้องปฏิบัติการมี IC50 ไม่สอดคล้อง | ปานกลาง | Noise ในตัวแปรตาม | Dedup ด้วย median + ตัดกลุ่มที่ SD > 1.0 log (Phase 2 ขั้นที่ 10) |
| **R6** | `pip install rdkit/shap` ล้มเหลวบน Python 3.13 | ปานกลาง | บล็อกทั้งโครงงาน | สร้าง venv ด้วย Python 3.11 เป็นแผนสำรอง — ทดสอบใน Phase 0 ทันที |
| **R7** | GridSearchCV ช้าเกินไปบน 4 cores | ปานกลาง | เสียเวลา | ใช้ `RandomizedSearchCV(n_iter=60)`, ตั้ง `n_jobs=4`, ลด grid ของ SVR ก่อนเพราะ SVR ช้าที่สุดกับ 1,028 features |
| **R8** | Fingerprint 1,024 bits ทำให้ SHAP ช้ามาก | ปานกลาง | Phase 7 ล่าช้า | ใช้ `TreeExplainer` (เร็ว) และ subsample test set 300–500 แถว |
| **R9** | ไม่ได้รับอนุมัติให้ดึงข้อมูลใหม่ | ต่ำ | ต้องเปลี่ยนแผน | ใช้แผน B ในหัวข้อ 4.1 (DDH 94 สาร + Repeated Nested CV) |
| **R10** | DDH ซ้ำกับ ChEMBL ทำให้ external validation ไม่ valid | ต่ำ | ผลลวง | ตรวจ overlap ด้วย InChIKey ก่อนรายงาน (T2.x) |
| **R11** | โฟลเดอร์ยังไม่เป็น git repository | ต่ำ | สูญเสียงาน | `git init` ใน Phase 0 |

---

## ส่วนที่ 11: ตารางเวลาดำเนินงาน

### 11.1 แผนโดยละเอียด (ปรับจากตารางที่ 1 ในหัวข้อ 9 ของ PDF)

| Phase | งาน | ระยะเวลา | ผลผลิตหลัก |
|---|---|---|---|
| 0 | เตรียมสภาพแวดล้อม + Data Audit + ขออนุมัติ | 2–3 วัน | `requirements.txt`, `tbl01` |
| 1 | Data Acquisition (ChEMBL + DDH) | 3–4 วัน | `chembl_*_raw.csv` |
| 2 | Data Cleaning & Standardization | 4–5 วัน | **`clean_bioactivity_master.csv`** |
| 3 | Feature Extraction (1,028 features) | 2–3 วัน | `features_hybrid.parquet`, `bitinfo.pkl` |
| 4 | EDA | 3 วัน | fig01–fig02c |
| 5 | Pipeline + 5-Fold CV (E1–E5) | 4–5 วัน | `tbl02_cv_results.csv` |
| 6 | GridSearchCV + Test Evaluation | 4 วัน | `tbl03_test_results.csv`, fig03–fig04b |
| 7 | SHAP + SAR Bit Mapping | 5–6 วัน | fig05–fig07, ตาราง SAR |
| 8 | External Validation + รายงาน + นำเสนอ | 5 วัน | `final_report.md`, สไลด์ |
| | **รวม** | **≈ 32–38 วันทำการ (7–8 สัปดาห์)** | |

### 11.2 การจับคู่กับ Gantt เดิม (มิ.ย.–พ.ย. 2569)

| งานใน Gantt เดิม | สถานะ | Phase ในแผนนี้ |
|---|---|---|
| 1. กำหนดหัวข้อโครงงาน | ✅ เสร็จ | — |
| 2. ศึกษาค้นคว้าทฤษฎีที่เกี่ยวข้อง | ✅ เสร็จ | — |
| 3. กำหนดขอบเขตและเป้าหมาย | ✅ เสร็จ | — |
| 4. หาชุดข้อมูล | 🔄 **ต้องทำใหม่** | Phase 0–1 |
| 5. จัดโครงงาน (ก.ย.–พ.ย.) | ⏳ กำลังเริ่ม | Phase 2–7 |
| 6. นำเสนอผลการทดลอง (พ.ย.) | ⏳ | Phase 8 |

> **จุดวิกฤต (Critical Path):** งานข้อ 4 "หาชุดข้อมูล" ตาม Gantt ควรเสร็จตั้งแต่ ส.ค.–ก.ย. แต่จากการตรวจสอบพบว่าข้อมูลที่มียังไม่ตรงวัตถุประสงค์ **จึงควรเร่ง Phase 0–1 ให้เสร็จภายในสัปดาห์แรกของเดือนกันยายน** เพื่อไม่ให้กระทบกำหนดนำเสนอเดือนพฤศจิกายน

### 11.3 Milestones

| M# | Milestone | กำหนด | เกณฑ์ผ่าน |
|---|---|---|---|
| M1 | ข้อมูลพร้อมใช้ | สิ้น Phase 2 | `clean_bioactivity_master.csv` ผ่าน validation ครบทุกข้อ |
| M2 | Baseline model ทำงาน | สิ้น Phase 5 | ทุกโมเดลชนะ Dummy อย่างมีนัยสำคัญ |
| M3 | ผลลัพธ์สุดท้าย freeze | สิ้น Phase 6 | Test set ถูกใช้แล้วครั้งเดียว, โมเดลบันทึกแล้ว |
| M4 | SAR อธิบายได้ | สิ้น Phase 7 | ≥ 6 บิตพร้อมภาพและคำอธิบายเชิงเคมี |
| M5 | ส่งมอบ | สิ้น Phase 8 | รายงาน + สไลด์ + reproducibility check ผ่าน |

---

## ภาคผนวก ก: คำสั่งและสคริปต์อ้างอิง

### ก.1 ติดตั้ง Dependencies

```bash
python -m venv .venv
.venv\Scripts\activate
pip install --upgrade pip
pip install pandas numpy scikit-learn matplotlib seaborn joblib tqdm pyarrow
pip install rdkit xgboost shap chembl-webresource-client
pip freeze > requirements.txt
```

### ก.2 ตรวจสอบ target ChEMBL (ยืนยันแล้วเมื่อ 8 ก.ย. 2569)

```python
import json, urllib.request
base = "https://www.ebi.ac.uk/chembl/api/data/"

# ยืนยัน target
t = json.load(urllib.request.urlopen(base + "target/CHEMBL4523582?format=json"))
print(t["pref_name"], "|", t["organism"], "|", t["target_type"])
# → Replicase polyprotein 1ab | Severe acute respiratory syndrome coronavirus 2 | SINGLE PROTEIN

# นับจำนวน activities
q = "activity?target_chembl_id=CHEMBL4523582&standard_type=IC50&format=json&limit=1"
print(json.load(urllib.request.urlopen(base + q))["page_meta"]["total_count"])   # 6220
print(json.load(urllib.request.urlopen(
    base + q + "&standard_units=nM"))["page_meta"]["total_count"])                # 6050
```

### ก.3 ดึงข้อมูลทั้งหมด (paging)

```python
import pandas as pd, json, urllib.request, time

def fetch_all(target_id="CHEMBL4523582", limit=1000):
    base = "https://www.ebi.ac.uk/chembl/api/data/activity"
    url = (f"{base}?target_chembl_id={target_id}&standard_type=IC50"
           f"&standard_units=nM&standard_relation==&format=json&limit={limit}")
    rows = []
    while url:
        d = json.load(urllib.request.urlopen(url, timeout=120))
        rows.extend(d["activities"])
        nxt = d["page_meta"]["next"]
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
        time.sleep(0.3)                      # เคารพ rate limit ของ EBI
    return pd.DataFrame(rows)

df = fetch_all()
df.to_csv("data/raw/chembl_CHEMBL4523582_ic50_raw.csv", index=False)
```

### ก.4 กรอง assay เฉพาะ 3CL-Pro (T1.4)

```python
mask = df["assay_description"].str.contains(
    r"3CL|3C-like|Mpro|M\s*pro|main\s+protease", case=False, na=False)
print(f"เก็บไว้ {mask.sum()} / {len(df)} records")
df_mpro = df[mask].copy()
```

### ก.5 Standardization และคำนวณ pIC50

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

lfc  = rdMolStandardize.LargestFragmentChooser()
norm = rdMolStandardize.Normalizer()
unch = rdMolStandardize.Uncharger()

def standardize(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    m = unch.uncharge(norm.normalize(lfc.choose(m)))
    return Chem.MolToSmiles(m, canonical=True), Chem.MolToInchiKey(m)

df_mpro["pIC50"] = 9 - np.log10(df_mpro["standard_value"].astype(float))
```

### ก.6 สกัด Hybrid Features 1,028 มิติ (T3.1–T3.5)

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdFingerprintGenerator

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

def featurize(smi):
    m  = Chem.MolFromSmiles(smi)
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitInfoMap()
    fp = gen.GetFingerprintAsNumPy(m, additionalOutput=ao)
    desc = [Descriptors.MolWt(m), Descriptors.MolLogP(m),
            Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m)]
    return fp, desc, ao.GetBitInfoMap()

FP_COLS   = [f"FP_{i}" for i in range(1024)]
DESC_COLS = ["MW", "LogP", "HBD", "HBA"]      # รวม 1,028 features
```

### ก.7 Pipeline ที่ป้องกัน Data Leakage (หัวใจของ O2)

```python
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler

def make_pipe(estimator):
    prep = ColumnTransformer([
        ("fp",   "passthrough",    FP_COLS),
        ("desc", StandardScaler(), DESC_COLS),
    ], remainder="drop")
    return Pipeline([("prep", prep), ("model", estimator)])
```

### ก.8 Scaffold Split (E2)

```python
from collections import defaultdict
from rdkit.Chem.Scaffolds import MurckoScaffold

def scaffold_split(smiles_list, test_frac=0.25):
    buckets = defaultdict(list)
    for i, s in enumerate(smiles_list):
        sc = MurckoScaffold.MurckoScaffoldSmiles(smiles=s, includeChirality=False)
        buckets[sc].append(i)
    groups = sorted(buckets.values(), key=len, reverse=True)   # กลุ่มใหญ่เข้า train ก่อน
    n_train = len(smiles_list) - int(len(smiles_list) * test_frac)
    train, test = [], []
    for g in groups:
        (train if len(train) + len(g) <= n_train else test).extend(g)
    return train, test
```

### ก.9 SHAP + Bit Mapping (O3)

```python
import shap, pickle, numpy as np
from rdkit import Chem
from rdkit.Chem import Draw

best  = pipe.named_steps["model"]
Xt    = pipe.named_steps["prep"].transform(X_test)
names = FP_COLS + DESC_COLS

explainer   = shap.TreeExplainer(best)
shap_values = explainer.shap_values(Xt)
mean_abs    = np.abs(shap_values).mean(axis=0)          # Importance_j ตามสูตร 4.5

top20 = np.argsort(mean_abs)[::-1][:20]
for j in top20:
    print(names[j], round(mean_abs[j], 4))

# วาดโครงสร้างย่อยของบิตสำคัญ
bitinfo = pickle.load(open("data/processed/bitinfo.pkl", "rb"))

def draw_bit(smiles, bit_id, out_png):
    m   = Chem.MolFromSmiles(smiles)
    img = Draw.DrawMorganBit(m, bit_id, bitinfo[smiles])
    img.save(out_png)
```

---

## เอกสารอ้างอิงเพิ่มเติมจากที่ระบุใน PDF

| # | อ้างอิง | ใช้สนับสนุนส่วนใด |
|---|---|---|
| A1 | ChEMBL Web Services — https://www.ebi.ac.uk/chembl/api/data/docs | Phase 1 (T1.2) |
| A2 | Bemis & Murcko (1996). *J. Med. Chem.* 39(15), 2887–2893 — Molecular frameworks | Scaffold Split (E2) |
| A3 | Sheridan (2013). *J. Chem. Inf. Model.* 53(4), 783–790 — Time-split / realistic validation | เหตุผลที่ต้องทำ E2 |
| A4 | Landrum, G. — RDKit Documentation (`DrawMorganBit`, `MolStandardize`) | Phase 2–3, 7 |
| A5 | Lundberg & Lee (2017). NIPS — A Unified Approach to Interpreting Model Predictions | Phase 7 |
| A6 | Maggiora (2006). *J. Chem. Inf. Model.* 46(4), 1535 — Activity cliffs | ความเสี่ยง R3 |

---

*เอกสารนี้จัดทำจากการอ่านข้อเสนอโครงงานฉบับเต็ม (10 หน้า) และการตรวจสอบชุดข้อมูลจริงทั้ง 5 ไฟล์ใน `Data_set/` เมื่อวันที่ 8 กันยายน 2569*

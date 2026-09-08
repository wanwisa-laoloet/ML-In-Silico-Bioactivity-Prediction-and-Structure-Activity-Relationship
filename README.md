# SilicoBioactivityPredict

การประยุกต์ใช้การเรียนรู้ของเครื่องเพื่อทำนายฤทธิ์ทางชีวภาพและวิเคราะห์ความสัมพันธ์ระหว่างโครงสร้างกับฤทธิ์ของสารยับยั้ง **SARS-CoV-2 Main Protease (Mpro / 3CLpro)** แบบ In Silico

รายวิชา CP352101 วิทยาการข้อมูลขั้นแนะนำ · วิทยาลัยการคอมพิวเตอร์ มหาวิทยาลัยขอนแก่น

---

## ผลลัพธ์โดยย่อ

| รายการ | ค่า |
|---|---|
| ชุดข้อมูล | ChEMBL `CHEMBL4523582` (SARS-CoV-2 3CL-Pro) — 2,198 โมเลกุลไม่ซ้ำ |
| Features | 1,028 (Morgan ECFP4 1,024 bits + Lipinski 4 ตัว) |
| โมเดลที่ดีที่สุด | **XGBoost** |
| Test R² (random split) | **0.7034** |
| Test R² (scaffold split) | 0.5127 |
| Spearman ρ | 0.8415 |

รายละเอียดทั้งหมดอยู่ใน [`reports/final_report.md`](reports/final_report.md)
แผนการดำเนินงานอยู่ใน [`Docs/Project_Plan.md`](Docs/Project_Plan.md)

---

## การติดตั้ง

โปรเจกต์นี้มี virtual environment เฉพาะอยู่ที่ `.venv/` (Python 3.13.5) ติดตั้งครบแล้ว

### เปิดใช้งาน venv

**PowerShell**
```powershell
.\.venv\Scripts\Activate.ps1
```

**Command Prompt**
```cmd
.venv\Scripts\activate.bat
```

**Git Bash**
```bash
source .venv/Scripts/activate
```

### สร้างใหม่จากศูนย์ (กรณีย้ายเครื่อง)

```bash
python -m venv .venv
```

```bash
./.venv/Scripts/python.exe -m pip install -r requirements.txt
```

### แพ็กเกจที่ติดตั้ง

| Package | Version | ใช้ทำอะไร |
|---|---|---|
| rdkit | 2026.3.6 | จัดการโครงสร้างเคมี, Morgan FP, Lipinski, DrawMorganBit |
| scikit-learn | 1.7.2 | Pipeline, ColumnTransformer, CV, RF, SVR, metrics |
| xgboost | 3.4.1 | XGBoost Regressor (โมเดลที่ดีที่สุด) |
| shap | 0.52.0 | TreeExplainer สำหรับวิเคราะห์ SAR |
| pandas | 2.3.2 | จัดการข้อมูลตาราง |
| numpy | 2.2.6 | คำนวณเชิงตัวเลข |
| scipy | 1.16.3 | Mann-Whitney U, Spearman ρ |
| matplotlib | 3.11.1 | สร้างกราฟ |
| seaborn | 0.13.2 | ธีมและ statistical plots |
| pyarrow | 25.0.1 | อ่าน/เขียนไฟล์ parquet |
| joblib | 1.5.2 | บันทึก/โหลดโมเดล |

---

## การรันซ้ำทั้งหมด

รันตามลำดับจากโฟลเดอร์ `src/` (ทุกขั้นตอนใช้ `RANDOM_STATE = 42`)

```bash
cd src && ../.venv/Scripts/python.exe acquire.py
```

```bash
cd src && ../.venv/Scripts/python.exe clean.py
```

```bash
cd src && ../.venv/Scripts/python.exe features.py
```

```bash
cd src && ../.venv/Scripts/python.exe eda.py
```

```bash
cd src && ../.venv/Scripts/python.exe model_cv.py
```

```bash
cd src && ../.venv/Scripts/python.exe tune_eval.py
```

```bash
cd src && ../.venv/Scripts/python.exe explain.py
```

```bash
cd src && ../.venv/Scripts/python.exe external_validation.py
```

```bash
cd src && ../.venv/Scripts/python.exe make_report.py
```

> `acquire.py` ต้องต่ออินเทอร์เน็ตเพื่อเรียก ChEMBL Web Services
> `tune_eval.py` ใช้เวลานานที่สุด (~10 นาทีบน CPU 4 cores)

---

## โครงสร้างโปรเจกต์

```
├── .venv/                      # virtual environment เฉพาะโปรเจกต์
├── Data_set/                   # ข้อมูลดิบต้นฉบับ (READ-ONLY)
├── Docs/
│   ├── Machine Learning Application ....pdf   # ข้อเสนอโครงงาน
│   ├── Project_Plan.md         # แผนดำเนินงาน 12 ส่วน
│   └── execution_log.md        # บันทึกความคืบหน้าและการรันสคริปต์
├── data/
│   ├── raw/                    # ผลดิบจาก ChEMBL API
│   ├── interim/                # หลังกรอง assay
│   └── processed/              # clean_bioactivity_master.csv, features, predictions
├── src/                        # โค้ดทั้งหมด 11 ไฟล์
├── models/                     # best_model_xgboost_*.joblib
├── reports/
│   ├── figures/                # fig01–fig09 (18 ไฟล์)
│   ├── tables/                 # tbl01–tbl08 (10 ไฟล์)
│   ├── latex/                  # เล่มรายงานวิทยานิพนธ์ฉบับเต็ม (report.pdf)
│   ├── Full_Report_5Chapters.md # รายงาน 5 บทฉบับ Markdown
│   └── final_report.md         # รายงานสรุปผลฉบับย่อ
└── requirements.txt
```

---

## ข้อจำกัดสำคัญ

- เป็นการศึกษา **In Silico** เท่านั้น ไม่มีการยืนยันในห้องปฏิบัติการ
- SHAP แสดง **ความสัมพันธ์เชิงสถิติ** ไม่ใช่กลไกการจับเชิงสาเหตุระดับ 3 มิติ
- **External validation บนชุด DDH ให้ R² = −1.72** เพราะ 84.4% ของสารอยู่นอก
  Applicability Domain — ห้ามใช้ค่า pIC50 ที่ทำนายเป็นตัวเลขสัมบูรณ์กับสารนอกโดเมน

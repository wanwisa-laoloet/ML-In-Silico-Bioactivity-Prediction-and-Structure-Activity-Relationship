# โครงสร้างโปรเจกต์รายงาน LaTeX (LaTeX Report Project Structure)

รายงานโครงงานคอมพิวเตอร์ CP352101  
**การประยุกต์ใช้การเรียนรู้ของเครื่องเพื่อทำนายฤทธิ์ทางชีวภาพและวิเคราะห์ความสัมพันธ์ระหว่างโครงสร้างกับฤทธิ์ของสารยับยั้ง SARS-CoV-2 Main Protease แบบ In Silico**

---

## โครงสร้างไดเรกทอรี (Directory Layout)

```
reports/latex/
├── report.tex                    # ไฟล์หลัก (Main Master Document)
├── .gitignore                    # ละเว้นไฟล์ build ชั่วคราว (.aux, .log, .toc ฯลฯ)
├── .latexmkrc                    # ตั้งค่า latexmk (XeLaTeX engine)
├── build.bat                     # สคริปต์คอมไพล์ XeLaTeX 2 รอบสำหรับ Windows
├── clean.bat                     # สคริปต์ลบไฟล์ build artifacts ชั่วคราว
├── README.md                     # เอกสารแนะนำการใช้งาน
│
├── preamble/                     # การตั้งค่าแพ็กเกจและฟอนต์
│   ├── preamble.tex              # Geometry, fontspec (ฟอนต์ไทย), จัดระยะ, หัว/ท้ายกระดาษ
│   └── commands.tex              # นิยามคำสั่งย่อ (\Mpro, \pICfifty), สี, กล่องเน้นข้อความ
│
├── frontmatter/                  # ส่วนหน้าของรายงาน
│   ├── titlepage.tex             # หน้าปกใน
│   ├── abstract_th.tex           # บทคัดย่อภาษาไทย
│   ├── abstract_en.tex           # Abstract ภาษาอังกฤษ
│   └── acknowledgments.tex       # กิตติกรรมประกาศ
│
├── chapters/                     # เนื้อหาหลักแยกตามบท
│   ├── 01_intro.tex              # บทที่ 1 บทนำ
│   ├── 02_theory.tex             # บทที่ 2 ทฤษฎีและงานวิจัยที่เกี่ยวข้อง
│   ├── 03_method.tex             # บทที่ 3 วิธีการดำเนินงาน
│   ├── 04_results.tex            # บทที่ 4 ผลการทดลองและอภิปรายผล
│   └── 05_conclusion.tex         # บทที่ 5 สรุปผลการทดลองและข้อเสนอแนะ
│
└── backmatter/                   # ส่วนท้ายของรายงาน
    ├── bibliography.tex          # รายการบรรณานุกรม
    ├── references.bib            # ฐานข้อมูลบรรณานุกรม BibTeX (.bib)
    └── appendices/               # ภาคผนวก
        ├── appendix_a_code.tex   # ภาคผนวก ก โครงสร้างไฟล์และโค้ดโปรแกรม
        └── appendix_b_files.tex  # ภาคผนวก ข รายการไฟล์ภาพและตารางผลลัพธ์
```

---

## ข้อกำหนดระบบและฟอนต์ (Requirements & Fonts)

- **TeX Engine:** XeLaTeX (แนะนำผ่าน MiKTeX หรือ TeX Live)
- **ฟอนต์ภาษาไทย:** ระบบจะเลือกฟอนต์อัตโนมัติตามลำดับที่มีในเครื่อง:
  1. `TH Sarabun New` (มาตรฐานราชการ/วิทยานิพนธ์)
  2. `Leelawadee UI` (ฟอนต์ระบบ Windows มาตรฐาน)
  3. `Angsana New`
  4. `Tahoma`
- **ฟอนต์ Monospace:** `Consolas` หรือ `Courier New`

---

## วิธีการคอมไพล์ (How to Compile)

### วิธีที่ 1: ดับเบิลคลิกหรือรันผ่าน `build.bat` (แนะนำสำหรับ Windows)
```cmd
build.bat
```

### วิธีที่ 2: ใช้คำสั่ง XeLaTeX โดยตรง
```bash
xelatex report.tex
xelatex report.tex
```
*(ต้องรัน 2 ครั้ง เพื่อให้สารบัญและหมายเลขอ้างอิงตรงกัน)*

### วิธีที่ 3: ใช้ `latexmk`
```bash
latexmk -xelatex report.tex
```

---

## การล้างไฟล์ชั่วคราว (Clean Build Artifacts)

```cmd
clean.bat
```
หรือใช้ `latexmk -c`

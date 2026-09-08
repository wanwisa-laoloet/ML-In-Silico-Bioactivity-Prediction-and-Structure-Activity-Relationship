@echo off
chcp 65001 > nul
echo Cleaning LaTeX auxiliary and build files...

del /q /f *.aux *.bbl *.bcf *.blg *.fdb_latexmk *.fls *.lof *.log *.lot *.out *.run.xml *.synctex.gz *.toc *.xdv 2>nul
del /s /q /f *.aux 2>nul

echo Clean complete!

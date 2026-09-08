@echo off
chcp 65001 > nul
echo ========================================================
echo  Compiling LaTeX Report (XeLaTeX)
echo ========================================================

echo [Pass 1/2] Running XeLaTeX...
xelatex -interaction=nonstopmode -halt-on-error report.tex
if errorlevel 1 goto error

echo [Pass 2/2] Running XeLaTeX (resolving references and TOC)...
xelatex -interaction=nonstopmode -halt-on-error report.tex
if errorlevel 1 goto error

echo ========================================================
echo  Compilation successful: report.pdf
echo ========================================================
goto end

:error
echo ========================================================
echo  Compilation failed! Please check report.log for errors.
echo ========================================================
exit /b 1

:end

@echo off
set paper=FDS_Config_Management_Plan

Title Building %paper%

git describe --abbrev=7 --long --dirty > gitinfo.txt
set /p gitrevision=<gitinfo.txt
echo \newcommand^{\gitrevision^}^{%gitrevision%^} > ..\Bibliography\gitrevision.tex

echo pass 1
pdflatex -interaction nonstopmode %paper% > %paper%.err
bibtex %paper% > %paper%.err
echo pass 2
pdflatex -interaction nonstopmode %paper% > %paper%.err
echo pass 3
pdflatex -interaction nonstopmode %paper% > %paper%.err
echo pass 4
pdflatex -interaction nonstopmode %paper% > %paper%.err

python ..\scripts\check_manuals.py --file %paper%.tex --outname %paper%_py.err --suppressconsole

find "! LaTeX Error:" %paper%.err
find "Fatal error" %paper%.err
find "Error:" %paper%.err

find "Error:" %paper%_py.err
find "Warning:" %paper%_py.err
find "Misspelt" %paper%_py.err

echo %paper% build complete
pause


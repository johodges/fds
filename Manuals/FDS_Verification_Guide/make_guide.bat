@echo off
set paper=FDS_Verification_Guide

Title Building %paper%

git describe --abbrev=7 --long --dirty > gitinfo.txt
set /p gitrevision=<gitinfo.txt
echo \newcommand^{\gitrevision^}^{%gitrevision%^} > ..\Bibliography\gitrevision.tex

set AUXUSER=..\FDS_User_Guide\FDS_User_Guide.aux
set PDFUSER=..\FDS_User_Guide\FDS_User_Guide.pdf

if not exist %AUXUSER% goto else1
  copy %AUXUSER%
  goto endif1
:else1
  echo ***warning: %AUXUSER% does not exist. Build the FDS Users
  echo "            guide before building the Verification guide
:endif1


echo pass 1
pdflatex -interaction nonstopmode %paper% > %paper%.err
bibtex %paper% > %paper%.err
echo pass 2
pdflatex -interaction nonstopmode %paper% > %paper%.err
echo pass 3
pdflatex -interaction nonstopmode %paper% > %paper%.err
echo pass 4
pdflatex -interaction nonstopmode %paper% > %paper%.err

if not exist %PDFUSER% goto endif2
  copy %PDFUSER%
:endif2

python ..\scripts\check_manuals.py --datafile ..\scripts\files_to_check_ver.txt --outname %paper%_py.err --suppressconsole

find "! LaTeX Error:" %paper%.err
find "Fatal error" %paper%.err
find "Error:" %paper%.err

find "Error:" %paper%_py.err
find "Warning:" %paper%_py.err
find "Misspelt" %paper%_py.err

echo %paper% build complete
pause


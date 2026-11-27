@echo off
setlocal EnableExtensions
title Dendra2 Analysis - Select workbook and run
cd /d "%~dp0"

echo ============================================================
echo  Dendra2 turnover analysis
echo ============================================================
echo.
echo  1. Select the Excel workbook.
echo  2. Wait while statistics and figures are generated.
echo  3. The report and results folder will open automatically.
echo.

set "RSCRIPT="
for /f "delims=" %%P in ('where Rscript.exe 2^>nul') do if not defined RSCRIPT set "RSCRIPT=%%P"
if not defined RSCRIPT for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RSCRIPT=%%D\bin\Rscript.exe"
if not defined RSCRIPT for /d %%D in ("%LOCALAPPDATA%\Programs\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RSCRIPT=%%D\bin\Rscript.exe"
if not defined RSCRIPT for /d %%D in ("%USERPROFILE%\R\R-*") do if exist "%%D\bin\x64\Rscript.exe" set "RSCRIPT=%%D\bin\x64\Rscript.exe"

if not defined RSCRIPT (
  echo R is not installed or could not be found.
  echo Download R for Windows from: https://cran.r-project.org/bin/windows/base/
  echo Install with the default options, then run this launcher again.
  echo.
  pause
  exit /b 1
)

set "WORKBOOK=%~1"
if not defined WORKBOOK (
  for /f "usebackq delims=" %%F in (`powershell.exe -NoProfile -STA -Command "Add-Type -AssemblyName System.Windows.Forms; $d=New-Object System.Windows.Forms.OpenFileDialog; $d.Title='Select the Dendra2 Excel workbook'; $d.Filter='Excel workbooks (*.xlsx)|*.xlsx'; $d.Multiselect=$false; if($d.ShowDialog() -eq [System.Windows.Forms.DialogResult]::OK){$d.FileName}"`) do set "WORKBOOK=%%F"
)

if not defined WORKBOOK (
  echo No workbook was selected. Nothing was changed.
  pause
  exit /b 1
)

echo Using R: %RSCRIPT%
echo Workbook: %WORKBOOK%
echo.
"%RSCRIPT%" --vanilla "%~dp0scripts\run_analysis.R" "%WORKBOOK%"
if errorlevel 1 goto failed

set "RESULTS_DIR="
if exist "%~dp0scripts\LAST_RESULTS_PATH.txt" set /p RESULTS_DIR=<"%~dp0scripts\LAST_RESULTS_PATH.txt"
echo.
echo ============================================================
echo  FINISHED SUCCESSFULLY
echo ============================================================
if defined RESULTS_DIR (
  echo Results: %RESULTS_DIR%
  start "" "%RESULTS_DIR%"
  if exist "%RESULTS_DIR%\Dendra2_Analysis_Report.html" start "" "%RESULTS_DIR%\Dendra2_Analysis_Report.html"
)
echo.
pause
exit /b 0

:failed
echo.
echo ============================================================
echo  ANALYSIS STOPPED
echo ============================================================
echo Read the final error above. Existing data were not modified.
echo Any completed diagnostic files remain in the timestamped results folder.
echo.
pause
exit /b 1

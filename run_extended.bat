@echo off
cd /d "%~dp0"
python reproduce.py --mode extended --out reproduced_results_extended
if errorlevel 1 exit /b 1

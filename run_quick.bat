@echo off
cd /d "%~dp0"
python reproduce.py --mode quick
if errorlevel 1 exit /b 1

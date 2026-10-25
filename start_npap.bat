@echo off
rem NPAP one-click launcher: starts the web server and opens the browser.
rem First run creates the venv and installs dependencies (a few minutes).
rem ASCII-only on purpose: cmd parses .bat files in the system ANSI codepage.
cd /d %~dp0

if exist .venv\Scripts\npap.exe goto run

echo [NPAP] First run: creating venv and installing dependencies (a few minutes)...
py -3.13 -m venv .venv
if errorlevel 1 python -m venv .venv
if errorlevel 1 (
    echo [NPAP] Python 3.13 not found. Install it from https://www.python.org/downloads/
    pause
    exit /b 1
)
.venv\Scripts\python -m pip install -e . -q

:run
echo [NPAP] Starting web server (http://127.0.0.1:8000 , Ctrl+C to stop)...
.venv\Scripts\npap.exe serve
rem Server exited (port occupied / crash / Ctrl+C): keep the window open so
rem the message above stays readable.
echo.
echo [NPAP] Server stopped. If it failed to start, read the message above.
pause

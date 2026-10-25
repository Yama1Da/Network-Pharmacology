@echo off
rem NPAP: stop the running web server (frees port 8000).
rem ASCII-only on purpose: cmd parses .bat files in the system ANSI codepage.
set PID=
for /f "tokens=5" %%p in ('netstat -ano ^| findstr ":8000.*LISTENING"') do set PID=%%p
if "%PID%"=="" (
    echo [NPAP] No server running on port 8000. Nothing to stop.
    pause
    exit /b 0
)
echo [NPAP] Stopping server (PID %PID%)...
taskkill /F /PID %PID% >nul 2>&1
echo [NPAP] Stopped. Start again with start_npap.bat
pause

@echo off
rem NPAP: stop the running web server (frees port 8000).
rem ASCII-only on purpose: cmd parses .bat files in the system ANSI codepage.
setlocal
set "FOUND="
for /f "tokens=5" %%p in ('netstat -ano ^| findstr ":8000.*LISTENING"') do (
    set "FOUND=1"
    echo [NPAP] Stopping server (PID %%p^)...
    taskkill /F /PID %%p >nul 2>&1
)
if not defined FOUND (
    echo [NPAP] No server running on port 8000. Nothing to stop.
    pause
    exit /b 0
)

rem Give it a moment to release the port, then verify. A survivor means
rem access was denied (process owned by another user/session) - say so
rem instead of a false "Stopped".
ping -n 3 127.0.0.1 >nul
netstat -ano | findstr ":8000.*LISTENING" >nul
if not errorlevel 1 (
    echo [NPAP] WARNING: something still listens on port 8000 ^(access denied?^).
    echo [NPAP] Check  : netstat -ano ^| findstr ":8000.*LISTENING"
    echo [NPAP] Then   : taskkill /F /PID ^<pid^>   ^(or close the app using it^)
    pause
    exit /b 1
)
echo [NPAP] Stopped. Start again with start_npap.bat
pause

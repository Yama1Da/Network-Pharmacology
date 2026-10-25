@echo off
rem NPAP one-click launcher: starts the web server and opens the browser.
rem First run creates the venv and installs dependencies (a few minutes).
rem ASCII-only on purpose: cmd parses .bat files in the system ANSI codepage.
cd /d %~dp0

rem Default install source: Tsinghua PyPI mirror (fast in China).
rem Remove " -i ..." below to use the official PyPI instead.
set "NPAP_PIP_INDEX=-i https://pypi.tuna.tsinghua.edu.cn/simple"

if exist .venv\Scripts\npap.exe goto check

echo [NPAP] First run: creating venv and installing dependencies (a few minutes)...
py -3.13 -m venv .venv
if errorlevel 1 python -m venv .venv
if errorlevel 1 (
    echo [NPAP] Python 3.13 not found. Install it from https://www.python.org/downloads/
    pause
    exit /b 1
)
echo [NPAP] Installing dependencies. Do not close this window.
echo [NPAP] One dot per second - a growing line means it is working.
call :pip_with_heartbeat
if errorlevel 1 goto pip_failed
goto check

:pip_failed
echo.
echo [NPAP] Dependency install FAILED - the server would crash on start.
echo [NPAP] Common causes: unstable network or antivirus blocking new exe files.
if exist .venv\pip_out.log type .venv\pip_out.log
echo [NPAP] Retry the install manually, then run this launcher again:
echo     .venv\Scripts\python -m pip install -e . %NPAP_PIP_INDEX%
pause
exit /b 1

:check
rem Health check before start: an interrupted install can leave npap.exe in place
rem while dependencies or the WY package mapping are missing. Test the real entry
rem point itself (python -c "import WY" would false-pass via the project directory),
rem then the serve-time imports. Repair instead of crashing on start.
.venv\Scripts\npap.exe --help >nul 2>&1
if errorlevel 1 goto repair
.venv\Scripts\python -c "import uvicorn, fastapi, pandas" >nul 2>&1
if not errorlevel 1 goto run
:repair
echo [NPAP] Dependencies incomplete - repairing. One dot per second.
call :pip_with_heartbeat
if errorlevel 1 goto pip_failed
.venv\Scripts\npap.exe --help >nul 2>&1
if errorlevel 1 goto pip_failed

:run
echo [NPAP] Starting web server (http://127.0.0.1:8000 , Ctrl+C to stop)...
.venv\Scripts\npap.exe serve
rem Server exited (port occupied / crash / Ctrl+C): keep the window open so
rem the message above stays readable.
echo.
echo [NPAP] Server stopped. If it failed to start, read the message above.
pause
exit /b 0

:pip_with_heartbeat
rem Runs pip quietly in a background cmd and prints one dot per second until
rem it finishes, so the window never looks frozen. 30-minute safety cap.
rem pip output goes to .venv\pip_out.log; result marker is .venv\pip_rc.tmp.
if exist .venv\pip_rc.tmp del .venv\pip_rc.tmp >nul 2>&1
start "" /b cmd /c ".venv\Scripts\python -m pip install -e . -q %NPAP_PIP_INDEX% >.venv\pip_out.log 2>&1 && echo OK>.venv\pip_rc.tmp || echo FAIL>.venv\pip_rc.tmp"
set /a HB_N=0
:hb_loop
if exist .venv\pip_rc.tmp goto hb_done
set /a HB_N+=1
if %HB_N% geq 1800 goto hb_done
<nul set /p "=."
ping -n 2 127.0.0.1 >nul
goto hb_loop
:hb_done
echo.
set "PIP_RC="
if exist .venv\pip_rc.tmp set /p PIP_RC=<.venv\pip_rc.tmp
if exist .venv\pip_rc.tmp del .venv\pip_rc.tmp >nul 2>&1
if "%PIP_RC%"=="OK" exit /b 0
exit /b 1

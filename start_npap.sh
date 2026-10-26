#!/usr/bin/env bash
# NPAP one-click launcher for Linux / macOS: starts the web server.
# First run creates the venv and installs dependencies (a few minutes).
# Usage: ./start_npap.sh   (or: bash start_npap.sh)
# Official PyPI instead of the Tsinghua mirror: NPAP_PIP_INDEX= ./start_npap.sh
set -u
cd "$(dirname "$0")" || { printf '[NPAP] Cannot enter the script directory.\n'; exit 1; }

# Default install source: Tsinghua PyPI mirror (fast in China).
# Override with an empty value to use the official PyPI instead.
NPAP_PIP_INDEX="${NPAP_PIP_INDEX:--i https://pypi.tuna.tsinghua.edu.cn/simple}"

VENV_PY=".venv/bin/python"
VENV_NPAP=".venv/bin/npap"

log() { printf '[NPAP] %s\n' "$*"; }

# Find a Python 3.13 interpreter (project requires >=3.13,<3.14).
find_py() {
    if command -v python3.13 >/dev/null 2>&1; then
        command -v python3.13
        return 0
    fi
    if command -v python3 >/dev/null 2>&1; then
        if python3 -c 'import sys; raise SystemExit(0 if sys.version_info[:2] == (3, 13) else 1)' >/dev/null 2>&1; then
            command -v python3
            return 0
        fi
    fi
    return 1
}

# Run pip quietly in the background and print one dot per second, so the
# window never looks frozen. 30-minute safety cap. Result in .venv/pip_rc.tmp.
pip_with_heartbeat() {
    rm -f .venv/pip_rc.tmp
    "$VENV_PY" -m pip install -e . -q $NPAP_PIP_INDEX >.venv/pip_out.log 2>&1 &
    local pid=$!
    printf '[NPAP] Installing dependencies. One dot per second - a growing line means it is working.\n'
    local n=0
    while [ ! -f .venv/pip_rc.tmp ] && [ "$n" -lt 1800 ]; do
        printf '.'
        sleep 1
        n=$((n + 1))
    done
    wait "$pid"
    local rc=$?
    rm -f .venv/pip_rc.tmp
    echo
    return "$rc"
}

pip_failed() {
    echo
    log "Dependency install FAILED - the server would crash on start."
    log "Common causes: unstable network, or a venv without pip."
    if [ -f .venv/pip_out.log ]; then tail -n 20 .venv/pip_out.log; fi
    log "Retry the install manually, then run this launcher again:"
    log "    $VENV_PY -m pip install -e . $NPAP_PIP_INDEX"
    read -rp "[NPAP] Press Enter to exit..." _
    exit 1
}

if [ ! -x "$VENV_NPAP" ]; then
    log "First run: creating venv and installing dependencies (a few minutes)..."
    if ! PY="$(find_py)"; then
        log "Python 3.13 not found. Install it from https://www.python.org/downloads/"
        log "(macOS: brew install python@3.13 | Debian/Ubuntu: apt install python3.13 python3.13-venv)"
        exit 1
    fi
    log "Using Python: $PY"
    "$PY" -m venv .venv || { log "venv creation failed."; exit 1; }
    "$VENV_PY" -m pip --version >/dev/null 2>&1 || {
        log "The venv has no pip. Debian/Ubuntu: sudo apt install python3.13-venv, then delete .venv and retry."
        exit 1
    }
    pip_with_heartbeat || pip_failed
fi

# Health check before start: an interrupted install can leave npap in place
# while dependencies or the WY package mapping are missing. Test the real
# entry point itself, then the serve-time imports. Repair instead of crashing.
if ! "$VENV_NPAP" --help >/dev/null 2>&1 || ! "$VENV_PY" -c "import uvicorn, fastapi, pandas" >/dev/null 2>&1; then
    log "Dependencies incomplete - repairing."
    pip_with_heartbeat || pip_failed
    "$VENV_NPAP" --help >/dev/null 2>&1 || pip_failed
fi

# Port pre-check: "address already in use" would surface as a server traceback.
# Catch it with a readable hint instead (bash /dev/tcp probe, no extra tools).
if (exec 3<>"/dev/tcp/127.0.0.1/8000") 2>/dev/null; then
    log "Port 8000 is already in use - the server may already be running."
    log "Open http://127.0.0.1:8000 first, or run ./stop_npap.sh and start again."
    exit 1
fi

log "Starting web server (http://127.0.0.1:8000 , Ctrl+C to stop)..."
"$VENV_NPAP" serve

echo
log "Server stopped. If it failed to start, read the message above."

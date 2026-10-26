#!/usr/bin/env bash
# NPAP: stop the running web server (frees port 8000). Linux / macOS.
set -u
PORT=8000

log() { printf '[NPAP] %s\n' "$*"; }

# Find PIDs listening on $PORT. Tools differ per platform:
# lsof is built into macOS and present on most Linux distros; fuser/ss are
# Linux fallbacks (ss cannot always show the pid without root).
find_listener_pids() {
    if command -v lsof >/dev/null 2>&1; then
        lsof -ti tcp:"$PORT" -s TCP:LISTEN 2>/dev/null
    elif command -v fuser >/dev/null 2>&1; then
        fuser -n tcp "$PORT" 2>/dev/null | tr -s ' \t' '\n\n' | grep -E '^[0-9]+$'
    elif command -v ss >/dev/null 2>&1; then
        ss -ltnp 2>/dev/null | awk -v p=":$PORT" '$4 ~ p { if (match($0, /pid=[0-9]+/)) print substr($0, RSTART+4, RLENGTH-4) }'
    fi
}

PIDS="$(find_listener_pids)"
if [ -z "$PIDS" ]; then
    log "No server running on port $PORT. Nothing to stop."
    exit 0
fi

log "Stopping server (PID $(printf '%s ' $PIDS))..."
kill $PIDS 2>/dev/null

# Give it 5 seconds to exit gracefully, then force-kill any survivor.
local_wait=0
while [ "$local_wait" -lt 5 ] && [ -n "$(find_listener_pids)" ]; do
    sleep 1
    local_wait=$((local_wait + 1))
done
PIDS="$(find_listener_pids)"
if [ -n "$PIDS" ]; then
    kill -9 $PIDS 2>/dev/null
    sleep 1
fi

# Verify the port is actually freed; a survivor means access was denied
# (process owned by another user) - say so instead of a false "Stopped".
if [ -n "$(find_listener_pids)" ]; then
    log "WARNING: something still listens on port $PORT (access denied? another user's process?)."
    log "Inspect with: lsof -i tcp:$PORT   then kill it manually (sudo if needed)."
    exit 1
fi

log "Stopped. Start again with ./start_npap.sh"

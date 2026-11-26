#!/usr/bin/env bash
# 8 MC-Läufe je Variante (kelvin_cap, tolman), 4 parallel; danach Mittelung mit runningAver.py
set -euo pipefail
cd "$(dirname "$0")"
N=8
PY=../venv/bin/python
for v in kelvin_cap tolman; do for i in $(seq 1 $N); do echo "$v $i"; done; done |
  xargs -P 4 -n 2 bash -c 'cd "$0" && ../'"$PY"' SRC/scripts/SPMC-0D-main-200326.py L500.case "$1" > run_$1.log 2>&1 && echo "fertig: $0 $1"'
for v in kelvin_cap tolman; do (cd $v && ../$PY SRC/scripts/runningAver.py L500.case $N > runningAver.log 2>&1 && echo "gemittelt: $v"); done

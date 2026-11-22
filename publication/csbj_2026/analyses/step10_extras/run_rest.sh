#!/bin/bash
cd "$(dirname "$0")"; PY=$HOME/igea-venv/bin/python
for d in HIV Hypoxia Vaccination; do $PY s10c_annotation_matched.py "$d" > "s10c_$d.log" 2>&1; done
$PY s10d_uncertainty.py > s10d.log 2>&1
$PY s10a_sensitivity.py > s10a.log 2>&1
echo finished > run_all.done

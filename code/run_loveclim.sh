#!/bin/bash
# LOVECLIM runs (needs data/external/ibs_hds from the IBS server; see data/external/IBS_HDS_README.md).
# Baseline: both climates, both chronologies, c_mix 0. Then the costed mixed strategy on LOVECLIM only.
set -e
cd "$(dirname "$0")/../work"
PY=../.venv/bin/python
N=${N:-400}
[ -f cross.npz ] || $PY ../code/build.py
$PY ../code/mc_ibs.py $N
for c in 0.1 0.2; do OOA_CLIMS=ibs OOA_C_MIX=$c OOA_TAG=_c$c $PY ../code/mc_ibs.py $N; done
$PY ../code/analyze_ibs.py
$PY ../code/review_compare_ibs.py

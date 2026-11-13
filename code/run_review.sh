#!/bin/bash
# Review reruns. Gate: c_mix=0, wait=0 must reproduce results/ before any sweep runs.
set -e
cd "$(dirname "$0")/../work"
PY=../.venv/bin/python
N=${N:-400}
$PY ../code/build.py
run() { [ -f mc_$1.pkl ] || $PY ../code/mc.py $N --c-mix $2 --wait $3 --tag $1; $PY ../code/exp.py $1 > /dev/null; }
run c0_w0 0 0
$PY ../code/review_compare.py gate      # exits non-zero if the baseline does not reproduce
for c in 0.05 0.1 0.2; do run c${c}_w0 $c 0; done
for w in 2 5 10; do run c0_w${w} 0 $w; done
run c0.1_w5 0.1 5
$PY ../code/exp.py c0_w0 --no-nefud > /dev/null
$PY ../code/exp.py c0.1_w5 --no-nefud > /dev/null
$PY ../code/review_compare.py summary

#!/bin/bash
# Production runs: each set at three mesh densities (m=1,2,3) for Lfac=4, and m=2,3 for Lfac=8, m=2 for Lfac=2.
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=<home>/miniconda3/envs/fenicsx-env/bin/python
cd "$(dirname "$0")"
for s in V1 S1 S2 S3 S4; do
  for ml in "1 4" "2 4" "3 4" "2 8" "3 8" "2 2" "4 4"; do
    set -- $ml
    $PY bilayer.py $s $1 $2 > logs/bilayer_${s}_m$1_L$2.txt 2>&1 || echo "FAILED $s $ml"
  done
done
echo ALLDONE

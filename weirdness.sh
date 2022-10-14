#!/bin/bash

cd build
make -j
rm results
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
echo "Variant,Component,Time,Size,L3 Misses,L3 Accesses,L3 Miss Ratio" > results

for FACTOR in 0.5 1 1.25 1.5 1.75 2 2.5; do
for EXEC in original hand right; do
  ./bin/2mm-$EXEC.exe $FACTOR >> results
done
done


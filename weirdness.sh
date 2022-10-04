#!/bin/bash

cd build
make -j
rm results
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
echo "Variant,Component,Time,Size" > results

for EXEC in original hand right; do
for FACTOR in 0.5 0.7 0.8 0.9 1 1.2 1.3 1.5 1.55 1.6 1.65 1.7 2 4; do
  ./bin/2mm-$EXEC.exe $FACTOR >> results
done
done


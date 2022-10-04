#!/bin/bash

cd build
make -j
rm results
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
echo "Variant,Component,Time,Size" > results

for EXEC in original hand; do
for FACTOR in 0.5 0.7 0.8 0.9 1 1.2 1.3 1.5 1.7; do
  ./bin/2mm-$EXEC.exe $FACTOR >> results
done
done


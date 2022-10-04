#!/bin/bash

cd build
make -j
rm results
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
echo "Variant,Component,Time,Size" > results

for EXEC in original hand; do
for FACTOR in 1 2 3 4 5 6; do
  ./bin/2mm-$EXEC.exe $FACTOR >> results
done
done


#!/bin/bash

mkdir build
cd build
cmake -DENABLE_OPENMP=On -DBLT_CXX_STD=c++17 ..
make -j8

objdump --disassemble bin/2mm-original.exe > original.ins
objdump --disassemble bin/2mm-hand.exe > hand.ins
objdump --disassemble bin/2mm-lock.exe > lock.ins


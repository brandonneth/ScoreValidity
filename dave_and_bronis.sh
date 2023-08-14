#!/bin/bash

# Script to run the model timing experiments

echo Running Linear Model Implementation
rm $HOME/libs/lassen/lib/libRAJA.a
rm -r $HOME/libs/lassen/include/RAJA

cd $HOME/RAJA
git checkout DissertationModelExperiment
cd lassen-build
rm -rf ./*
module load gcc/8.3.1 cmake/3.23.1
cmake ..  -DENABLE_OPENMP=On -DBLT_CXX_STD=c++17 -DCMAKE_INSTALL_PREFIX=$HOME/libs/lassen -DENABLE_TESTS=Off -DENABLE_EXAMPLES=Off -DRAJA_ENABLE_EXERCISES=Off -DCMAKE_BUILD_TYPE=RELWITHDEBINFO
make
if [[ $? -ne 0 ]]; then
	echo "RAJA Build failed."
	exit
fi

make install
if [[ $? -ne 0 ]]; then
	echo "RAJA install failed."
	exit
fi

cd $HOME/ScoreValidity
cd build
rm -rf ./*
cmake ..  -DENABLE_OPENMP=On -DBLT_CXX_STD=c++17 -DCMAKE_INSTALL_PREFIX=$HOME/libs/lassen -DENABLE_TESTS=Off -DENABLE_EXAMPLES=Off -DRAJA_ENABLE_EXERCISES=Off -DCMAKE_BUILD_TYPE=RELWITHDEBINFO
make -j
if [[ $? -ne 0 ]]; then
	echo "ScoreValidity build failed."
	exit
fi

cd $HOME/ScoreValidity
build/bin/dave_and_bronis.exe > LinearModelDB2.csv




echo Running optimized model

rm $HOME/libs/lassen/lib/libRAJA.a
rm -r $HOME/libs/lassen/include/RAJA

cd $HOME/RAJA
git checkout formatdecisions
cd lassen-build
rm -rf ./*
module load gcc/8.3.1 cmake/3.23.1
cmake ..  -DENABLE_OPENMP=On -DBLT_CXX_STD=c++17 -DCMAKE_INSTALL_PREFIX=$HOME/libs/lassen -DENABLE_TESTS=Off -DENABLE_EXAMPLES=Off -DRAJA_ENABLE_EXERCISES=Off -DCMAKE_BUILD_TYPE=RELWITHDEBINFO
make
make install

cd $HOME/ScoreValidity
cd build
rm -rf ./*
cmake ..  -DENABLE_OPENMP=On -DBLT_CXX_STD=c++17 -DCMAKE_INSTALL_PREFIX=$HOME/libs/lassen -DENABLE_TESTS=Off -DENABLE_EXAMPLES=Off -DRAJA_ENABLE_EXERCISES=Off -DCMAKE_BUILD_TYPE=RELWITHDEBINFO
make -j

cd $HOME/ScoreValidity
build/bin/dave_and_bronis.exe > NonLinearModelDB2.csv


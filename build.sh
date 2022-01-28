
module load clang/13.0.0
export CXX=clang++

mkdir build
cd build
rm CMakeCache.txt
cmake -DENABLE_OPENMP=On -DBLT_CXX_STD=c++17 -DENABLE_TESTS=Off -DENABLE_EXAMPLES=Off ..
make -j8

cd ..


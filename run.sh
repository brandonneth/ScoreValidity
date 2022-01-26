
module load clang/13.0.0
export CXX=clang++

mkdir build
cd build
rm CMakeCache.txt
cmake -DENABLE_OPENMP=On -DBLT_CXX_STD=c++17 -DENABLE_TESTS=Off -DENABLE_EXAMPLES=Off ..
make -j8

cd ..

./build/bin/benchmark1.exe > benchmark1.times
./build/bin/benchmark2.exe > benchmark2.times
./build/bin/benchmark3.exe > benchmark3.times
./build/bin/intro-example.exe > intro-example.times

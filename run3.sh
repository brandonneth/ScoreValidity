
module load clang/13.0.0
export CXX=clang++

./build.sh

echo "Running Benchmark 3"
./build/bin/benchmark3.exe > benchmark3.times


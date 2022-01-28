
module load clang/13.0.0
export CXX=clang++

./build.sh

echo "Running Benchmark 1"
./build/bin/benchmark1.exe > benchmark1.times
echo "Running Benchmark 2"
./build/bin/benchmark2.exe > benchmark2.times
echo "Running Benchmark 3"
./build/bin/benchmark3.exe > benchmark3.times


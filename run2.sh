
module load clang/13.0.0
export CXX=clang++

./build.sh

echo "Running Intro Example"
./build/bin/intro-example.exe > intro-example.times
echo "Running 2MM Experiment"
./build/bin/2mm-all.exe > 2mm-all-choices.txt

REM uses Intel compiler (ifort) on Windows.
REM assumes you already used Intel compilers to build MUMPS.

del build_intel\CMakeCache.txt

set FC=ifort
set CC=icl
set CXX=icpc

cmake -B build_intel -DCMAKE_BUILD_TYPE=Debug -DMUMPS_ROOT=c:/lib_intel/mumps-5.2.1 -DPython3_ROOT=c:/miniconda3

cmake --build build_intel --parallel

cd build_intel

ctest -V
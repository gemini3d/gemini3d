REM uses Intel compiler (ifort) on Windows.
REM assumes you already used Intel compilers to build MUMPS.

del build_intel\CMakeCache.txt

cmake -B build_intel -DCMAKE_BUILD_TYPE=Debug -DMUMPS_ROOT=c:/lib_intel/mumps-5.2.1

cmake --build build_intel --parallel

cd build_intel

ctest -V
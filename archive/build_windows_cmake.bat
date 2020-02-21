REM uses Intel compiler (ifort) on Windows.
REM assumes you already used Intel compilers to build MUMPS.
REM timeout statements are due to glitches with CMake on Windows and parallel builds causing truncated commands

del build_intel\CMakeCache.txt

set FC=ifort
set CC=icl
set CXX=icl

cmake -B build_intel -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug -DMUMPS_ROOT=c:/lib_intel/mumps-5.2.1 -DPython3_ROOT=c:/miniconda3
if %errorlevel% neq 0 exit /b %errorlevel%

cmake --build build_intel --parallel
if %errorlevel% neq 0 exit /b %errorlevel%

timeout 1
cd build_intel

ctest -V
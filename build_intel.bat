REM for Intel compiler on Windows
@echo off

del %~dp0\objects\CMakeCache.txt

cmake -G "MinGW Makefiles" -DMUMPS_ROOT=C:\mumps-5.1.2 ..

cmake --build -B %~dp0\objects -S %~dp0 -j

cd %~dp0\objects

ctest --output-on-failure

cd %~dp0

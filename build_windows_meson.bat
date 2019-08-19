REM uses Intel compiler (ifort) on Windows.

set FC=ifort
set CC=icl
set CXX=icpc

del build_intel\build.ninja

meson setup build_intel --buildtype release -DMUMPS_ROOT=c:/lib_intel/mumps-5.2.1

meson test -C build_intel -v
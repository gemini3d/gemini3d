REM uses Intel compiler (ifort) on Windows.
REM Using MPI on Windows requires special permissions by allowing MPI to login itself

set FC=ifort
set CC=icl
set CXX=icpc

del build_intel\build.ninja

meson setup build_intel -Dmkl_root=%MKLROOT%

meson test -C build_intel -v
@echo off

REM uses Intel compiler (ifort) on Windows.
REM assumes you already used Intel compilers to build MUMPS.


REM %userprofile%\intel.bat

python build.py intel --arg="-DMUMPS_ROOT=d:/mumps-5.1.2"

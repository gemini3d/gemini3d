REM for Intel compiler on Windows
@echo off

python build.py intel --args="-DMUMPS_ROOT=C:\mumps-5.1.2" -test

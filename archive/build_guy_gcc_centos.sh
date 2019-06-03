#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

set -e
set -u

PREFIX=/home/ggrubbs1/source

#======================================================
MPIPREFIX=
LAPACKPREFIX=$PREFIX/lapack-3.8.0
SCALAPACKPREFIX=$PREFIX/scalapack-2.0.2
MUMPSPREFIX=$PREFIX/MUMPS_5.1.1

export FC=~/.local/bin/mpifort
export CC=~/.local/bin/mpicc

# ============================================================
# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
LD_LIBRARY_PATH=

. script_utils/check.sh

. script_utils/build.sh


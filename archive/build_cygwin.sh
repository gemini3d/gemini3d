#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

PREFIX=/cygdrive/c/Users/ggrubbs1/Desktop/lib

#======================================================
MPIPREFIX=
LAPACKPREFIX=
SCALAPACKPREFIX=
MUMPSPREFIX=$PREFIX/MUMPS_5.1.1


export FC=mpifort
export CC=mpicc

# ============================================================
# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
LD_LIBRARY_PATH=

# ----- build ! -------
. script_utils/check.sh

. script_utils/build.sh

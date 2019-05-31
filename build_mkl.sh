#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

PREFIX=$HOME/lib_gemini_gcc8mkl

#======================================================
MPIPREFIX=
LAPACKPREFIX=
SCALAPACKPREFIX=
MUMPSPREFIX=$PREFIX/mumps-$SUFFIX

. script_utils/check.sh

. script_utils/build.sh

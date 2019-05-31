#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for ifort
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

PREFIX=$HOME/lib_gemini_intel18

FC=mpifort
CC=mpicc

. script_utils/ifort_self.sh

. script_utils/check.sh

. script_utils/build.sh


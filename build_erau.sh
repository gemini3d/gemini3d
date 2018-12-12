#!/bin/bash
#
# "-d" option makes this a Debug build
#
# This is for Matt's machine, so he can use the older libraries he's been working with for some time.

OPTS="-DMUMPS_ROOT=~/lib/MUMPS_4.10.0 -DSCALAPACK_ROOT=~/lib/scalapack-2.0.2"

[[ $1 == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"

cmake $OPTS -B objects .

cmake --build objects -j

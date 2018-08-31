#!/bin/sh
#
# This is for Matt's machine, so he can use the older libraries he's been working with for some time.

(
cd objects

cmake -DMUMPS_ROOT=~/lib/MUMPS_4.10.0 \
      -DLAPACK95_ROOT=~/lib/LAPACK95 \
      -DSCALAPACK_ROOT=~/lib/scalapack-2.0.2 \
      ..
      
make
)

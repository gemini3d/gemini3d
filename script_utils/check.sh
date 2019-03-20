#!/bin/bash

# this is helper script, not meant to be run directly

set -u
set -e

for d in $MPIPREFIX $LAPACKPREFIX $SCALAPACKPREFIX $MUMPSPREFIX
do
  [[ -z $d ]] && continue
  [[ -d $d ]] || { echo "ERROR: $d not found"; exit 1; }
done

OPTS="-DMPI_ROOT=$MPIPREFIX ${OPTS:-}"
OPTS="-DLAPACK_ROOT=$LAPACKPREFIX $OPTS"
OPTS="-DSCALAPACK_ROOT=$SCALAPACKPREFIX $OPTS"
OPTS="-DMUMPS_ROOT=$MUMPSPREFIX $OPTS"
OPTS="-DUSEGLOW=true $OPTS"

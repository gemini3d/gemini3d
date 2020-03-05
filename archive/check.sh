#!/bin/bash

# this is helper script, not meant to be run directly

set -u
set -e

# We omit $MPIPREFIX because many HPC have it installed for the compiler already
# ditto for $LAPACKPREFIX
for d in $SCALAPACKPREFIX $MUMPSPREFIX
do
  [[ -z $d ]] && continue
  [[ -d $d ]] || { echo "ERROR: $d not found"; exit 1; }
done

OPTS="-DMPI_ROOT=$MPIPREFIX ${OPTS:-}"
OPTS="-DLAPACK_ROOT=$LAPACKPREFIX $OPTS"
OPTS="-DSCALAPACK_ROOT=$SCALAPACKPREFIX $OPTS"
OPTS="-DMUMPS_ROOT=$MUMPSPREFIX $OPTS"

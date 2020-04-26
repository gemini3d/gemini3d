#!/bin/bash -l

#$ -l h_rt=24:00:00

# this is total per-node RAM
#$ -l mem_per_core=3G

# MPI
#$ -pe omp 28

echo "NSLOTS = $NSLOTS"

. $HOME/gcccompilers.sh

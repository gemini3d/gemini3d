#!/bin/bash -l

#$ -l h_rt=12:00:00

# this is total per-node RAM
#$ -l mem_per_core=3G

# MPI
#$ -pe mpi_8_tasks_per_node 16

echo "NSLOTS = $NSLOTS"

. $HOME/gcccompilers.sh

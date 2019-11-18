#!/bin/sh
# This is for Ubuntu only, it uses system libraries
# apt install libblacs-openmpi-dev libtmglib-dev liblapack-dev libblas-dev libmumps-dev

make -j LIBBLACS=-lblacs-openmpi TMG=-ltmglib LAPACK77=-llapack BLAS=-lblas MUMPS=-ldmumps -lmumps-common -lpord MUMPSDIR=/usr

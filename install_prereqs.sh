#!/bin/bash
# installs libraries
# For CentOS HPC systems, consider centos*.sh scripts that use "modules" or "extensions"

set -e 

case $OSTYPE in
linux*) 
  
  if [[ -f /etc/redhat-release ]]; then
     echo "This script is made for personal laptops/desktops. HPCs using CentOS have system-specific library setup."
     echo "If using gfortran, version >= 6 is required. Try devtoolset-7 if gcc/gfortran is too old on your system."
     yum install epel-release
     yum install pkg-config gcc-gfortran g++ cmake make 
     yum install MUMPS-openmpi-devel lapack-devel scalapack-openmpi-devel blacs-openmpi-devel openmpi-devel libmetis-devel libscotch-devel ptscotch-openmpi-devel atlas-devel
     yum install octave
  else
    apt update
    apt install gfortran g++ make  # repo cmake is often too old
    apt install libssl-dev unzip  # needed BEFORE building CMake
    apt install libmumps-dev liblapack-dev libscalapack-mpi-dev libblacs-mpi-dev libopenmpi-dev libmetis-dev libscotch-dev libptscotch-dev libatlas-base-dev 
    apt install --no-install-recommends octave
  fi
  ;;
darwin*)
  brew install cmake gcc make lapack scalapack openmpi octave
# MUMPS is a little more involved, but OK.
  brew tap dpo/openblas;
  brew tap-pin dpo/openblas;
  brew options mumps;
  brew install mumps;
  ;;
cygwin*)
  echo "please install prereqs via Cygwin setup.exe"
  echo "gcc-fortran make liblapack-devel libopenmpi-devel octave"
  echo "then use https://github.com/scivision/fortran-libs to get all other prereqs"
  ;;
esac


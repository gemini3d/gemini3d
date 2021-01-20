# Fortran 2008 compilers

GEMINI requires a Fortran compiler that handles standard Fortran syntax including

* submodule
* block

These compilers work easily with GEMINI3D on Linux, MacOS and Windows:

* Gfortran &ge; 7
* Intel oneAPI core + [HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html)

## Linux


```sh
apt install gfortran
# or
yum install gcc-gfortran
```

Typical HPC have the ability to switch to a recent GCC version without sudo.

## MacOS

Install the latest Gfortran via [Homebrew](https://brew.sh):

```sh
brew install gcc
```

Install everything needed for Gemini like:

```sh
brew install cmake gcc hdf5 openmpi lapack scalapack
```

### Ensure GCC is used instead of Clang

If using Homebrew, be sure Homebrew's GCC is used instead of AppleClang or other non-Homebrew compilers so that the Homebrew library ABIs match the compiler ABI.
By specifying a GCC version, actual GCC will be used instead of Clang.

```sh
FC=gfortran-10 CC=gcc-10 cmake -B build
```

If you need to specify MPI compiler wrappers, do like:

```sh
cmake -B build -DMPI_ROOT=~/lib_gcc/openmpi
```

## Windows

[MSYS2](https://www.scivision.dev/install-msys2-windows)
provides a comprehensive Windows development solution.
From the MSYS2 terminal, install GCC and Gfortran by:

```sh
pacman -S mingw-w64-x86_64-gcc-fortran
```

other Gemini3D-required libraries are obtained by:

```sh
pacman -S mingw-w64-x86_64-msmpi
pacman -S mingw-w64-x86_64-hdf5
pacman -S mingw-w64-x86_64-lapack
pacman -S mingw-w64-x86_64-scalapack
```

Install
[Microsoft MS-MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi-release-notes),
which gives `mpiexec`.

# Gemini3D compilers

GEMINI requires a Fortran compiler that handles standard Fortran syntax including "submodule" and "block".
GEMINI requires a C++17 compiler that handles [filesystem](https://en.cppreference.com/w/cpp/filesystem) stdlib.

These compilers are known to work with GEMINI3D on Linux, MacOS, and Windows:

* Gfortran (GCC): 7.5, 8.5, 9.3, 10.3, 11.1, 11.2
* Intel oneAPI 2021.x, 2022.x core + [HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html)
* Cray with GCC backend

Intel Parallel Studio XE (PSXE) 2020 (Intel Fortran 19.1) is replaced by no-cost oneAPI.

Some older point releases of GCC are known to be broken (example: GCC 7.4 and 8.1 are broken in general).
Currently, GCC 9 is the oldest version maintained by the GCC devs.

## Linux

HPC users usually can switch to a recent GCC version.
[GCC](./Linux_gcc.md) is an easy choice for Linux users.

[Intel oneAPI](./Linux_intel_oneapi.md)
provides Intel MPI, LAPACK, and Scalapack for Linux.

## MacOS

The Clang C and C++ compilers work fine with [MacOS GCC](./MacOS_gcc.md) in general.

[Intel oneAPI on MacOS](./MacOS_intel_oneapi.md) does not include MPI, so this repo will build MPI for oneAPI.

## Windows

Windows users can choose between
[Intel oneAPI](./Windows_intel_oneapi.md),
[MSYS2](./Windows_gcc.md),
[Windows Subsystem for Linux](./Linux_gcc.md), or Cygwin.

# Gemini3D compilers

GEMINI requires a Fortran compiler that handles standard Fortran syntax including "submodule" and "block".
GEMINI requires a C++ compiler that handles
[filesystem](https://en.cppreference.com/w/cpp/filesystem)
stdlib.

These compilers are known to work with GEMINI3D on Linux, MacOS, and Windows:

* Gfortran (GCC): 9.3, 9.5, 10.3, 11.1, 11.2, 11.3, 12.1, 12.2
* Intel oneAPI &ge; 2022.2, 2023.0 core + [HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html)
* Cray with GCC or Intel backend

Intel Parallel Studio XE is obsolete and replaced by no-cost Intel oneAPI.

Some older point releases of GCC are known to be broken (example: GCC 9.2 is broken in general).

## Linux

HPC users usually can switch to a recent GCC version.
[GCC](./Linux_gcc.md)
is an easy choice for Linux users.

[Intel oneAPI](./Linux_intel_oneapi.md)
is another choice that may yield better runtime speed with Intel CPUs.

## MacOS

The GCC C, C++, and Fortran compilers work fine.
Although the Clang, Clang++ compilers also work with Gfortran, if you experience segmentation faults on the self-tests, try GCC, G++ instead.

* [GCC](./MacOS_gcc.md): Apple Silicon or Intel CPU Macs
* [Intel oneAPI](./MacOS_intel_oneapi.md) for Intel CPU Macs **only**. Recommend GCC instead for most users.

## Windows

Windows users can choose between
[Windows Subsystem for Linux](./Linux_gcc.md),
[Intel oneAPI](./Windows_intel_oneapi.md),
or Cygwin.
Cygwin has noticeably slower performance in general for any program.

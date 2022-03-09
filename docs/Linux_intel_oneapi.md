# Build Gemini3D with Intel oneAPI on Linux

Intel oneAPI provides Intel MPI, LAPACK, and Scalapack on Linux.

## setup Intel oneAPI

Like LLVM, oneAPI relies on the underlying GCC compiler for libc and libstdc++.
Having too old GCC/libc will fail to build.
Ensure your system has GCC 8 or newer to work with Intel oneAPI on Linux.

Install
[oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)
with these options:

* Math Kernel Library (oneMKL)

Install
[oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)
with these options:

* Intel MPI library
* Intel C++ compiler
* Intel Fortran compiler

We recommend making a little shell script named like "~/intel_oneapi.sh".
The contents of this script would be like:

```sh
source /opt/intel/oneapi/setvars.sh

export FC=ifort CC=icx CXX=icpx
```

Intel LLVM "icx" and "icpx" replace the legacy "icc" and "icpc" compilers.

To enable oneAPI in this Terminal:

```sh
source ~/intel_oneapi.sh
```

## One-time Gemini3D external library setup

```sh
git clone https://github.com/gemini3d/external

cmake -S external -B external/build -DCMAKE_INSTALL_PREFIX=~/lib_intel

cmake --build external/build
```

that installs Gemini3d external libraries under ~/lib_intel.
This path is arbitrary but should be distinct between compilers.

## Build and Test Gemini3D

```sh
git clone https://github.com/gemini3d/gemini3d

cmake -S gemini3d -B build/gemini3d -G Ninja -DCMAKE_PREFIX_PATH=~/lib_intel

cmake --build gemini3d/build

ctest --test-dir gemini3d/build
```

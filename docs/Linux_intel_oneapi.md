# Build Gemini3D with Intel oneAPI on Linux

Intel oneAPI (no cost) provides Intel MPI, LAPACK, and Scalapack on Linux.

Like LLVM, oneAPI relies on the underlying GCC compiler for libc and libstdc++.
Having too old GCC/libc/libstdc++ will fail to build.
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

We recommend making a little shell script named like "~/oneapi.sh".
The contents of this script would be like:

```sh
source /opt/intel/oneapi/setvars.sh

export FC=ifx CC=icx CXX=icpx
```

"icx", "icpx", "ifx" are the LLVM-based Intel compilers for C, C++, and Fortran, respectively.

To enable oneAPI in this Terminal:

```sh
source ~/oneapi.sh

export CC=$CMPLR_ROOT/bin/icx
export CXX=$CMPLR_ROOT/bin/icpx
export FC=$CMPLR_ROOT/bin/ifx
```

GCC is used as a backend for oneAPI on Linux, and a too-old GCC won't work with C++ programs.
Use GCC 10 or newer to work with Intel oneAPI on Linux.

If there are problems with defaulting to old GCC, specify the GCC toolchain in "~/oneapi.sh" like:

```sh
export CXXFLAGS=--gcc-toolchain=/opt/rh/gcc-toolset-11/root/usr/
```

which can be determined like:

```sh
scl enable gcc-toolset-11 "which g++"
```

## One-time Gemini3D external library setup

```sh
git clone https://github.com/gemini3d/external

cmake -Dprefix=~/libgem_intel -P external/build.cmake
```

that installs Gemini3d external libraries under ~/libgem_intel.
This path is arbitrary but should be distinct between compilers.

## Build and Test Gemini3D

```sh
git clone https://github.com/gemini3d/gemini3d

cmake -S gemini3d -B build/gemini3d -DCMAKE_PREFIX_PATH=~/libgem_intel

cmake --build gemini3d/build

ctest --test-dir gemini3d/build
```

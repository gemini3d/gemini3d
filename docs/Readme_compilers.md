# Gemini3D compilers

GEMINI requires a Fortran compiler that handles standard Fortran syntax including "submodule" and "block".
GEMINI requires a C++17 compiler that handles [filesystem](https://en.cppreference.com/w/cpp/filesystem) stdlib.

These compilers are known to work with GEMINI3D on Linux, MacOS, and Windows:

* Gfortran (GCC): 7.5, 8.5, 9.3, 10.3, 11.1, 11.2
* Intel oneAPI 2021.x, 2022.x core + [HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html)

Intel Parallel Studio XE (PSXE) 2020 (Intel Fortran 19.1) is replaced by no-cost oneAPI.

Note: some older point releases of GCC are known to be broken (example: GCC 7.4 and 8.1 are broken in general).
Currently, GCC 9 is the oldest version maintained by the GCC devs.

## Linux

```sh
apt install g++ gfortran
# or
yum install gcc-c++ gcc-gfortran
```

HPC users usually can switch to a recent GCC version.

[Intel oneAPI](https://www.scivision.dev/intel-oneapi-fortran-install/)
provides Intel MPI, LAPACK, and Scalapack for Linux.

## MacOS

Install the latest Gfortran via [Homebrew](https://brew.sh):

```sh
brew install gcc
```

The Clang C and C++ compilers work fine with Gemini3D and Gfortran in general.

Intel oneAPI on MacOS does not include MPI, so you would have to build MPICH yourself.

## Windows

Windows users can choose between Intel oneAPI, MSYS2, Windows Subsystem for Linux, or Cygwin.

[MSYS2](https://www.scivision.dev/install-msys2-windows)
provides a comprehensive Windows development solution from the Windows terminal.
From the MSYS2 terminal, install GCC C++ and Fortran compilers:

```sh
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran
```

[Intel oneAPI](https://www.scivision.dev/intel-oneapi-fortran-install/)
provides Intel MPI, LAPACK, and Scalapack on Windows.
We do not use MSYS2/GCC libraries with Windows oneAPI as they are ABI incompatible.
Use the oneAPI Command Prompt on Windows.

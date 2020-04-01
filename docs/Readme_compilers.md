# Fortran 2008 compilers

Modern Fortran compilers are widely and freely available.
If you don't already have a Fortran compiler, it's straightforward to install.
GEMINI requires a Fortran compiler that handles Fortran 2008 syntax including

* submodule
* block

Most currently maintained Fortran compilers can do this.
Gfortran and Intel Fortran are two compilers known to work easily with GEMINI.
Here are details on a few Fortran compilers suitable for Gemini.

## Gfortran

Gfortran &ge; 6 has adequate Fortran 2008 support for Gemini.
Here are a few ways to get Gfortran on your computer:

### MacOS Homebrew

[Homebrew](https://brew.sh)
is a very popular development repository for MacOS.
Installing the latest GCC including Gfortran is simply:

```sh
brew install gcc
```

Install everything needed for Gemini like:

```sh
brew install cmake gcc openmpi lapack scalapack
```

If using Homebrew, be sure Homebrew's GCC is used instead of AppleClang or other non-Homebrew compilers so that the Homebrew library ABIs match the compiler ABI.

```sh
FC=gfortran-9 CC=gcc-9 cmake -B build
```

If you need to specify MPI compiler wrappers, do like:

```sh
cmake -B build -DMPI_ROOT=~/lib_gcc/openmpi-3.1.4
```

### Linux

Currently supported versions of Debian, Ubuntu, CentOS, and Windows Subsystem for Linux generally have Gfortran &ge; 6 available.
For CentOS, typical HPC will have the ability to switch to a recent GCC version--ask your HPC IT staff.

* Ubuntu / Debian / Windows Subsystem for Linux: `apt install gfortran`
* CentOS: `yum install gcc-gfortran`

### Windows

In general for Fortran development on Windows,
[MSYS2](https://www.scivision.dev/install-msys2-windows/)
provides a comprehensive development solution.
From the MSYS2 terminal, install GCC and Gfortran by:

```sh
pacman -S mingw-w64-x86_64-gcc-fortran
```

## Intel Fortran

Intel compilers are available at no-cost for academic instruction and open-source projects.
You need to have
[Intel Parallel Studio XE](https://software.intel.com/en-us/parallel-studio-xe)
that includes IntelMPI, BLACS and SCALAPACK.
Intel Fortran
[currently supported versions](https://software.intel.com/en-us/articles/intel-parallel-studio-xe-supported-and-unsupported-product-versions)
are generally targeted for Gemini support.

We regularly use the latest Intel compilers on Linux and Windows.
We haven't tried Intel compiler on MacOS, please let us know if it works or has issues.

### Intel MKL with Gfortran

Intel MKL with Gfortran should work, let us know if there's an issue.

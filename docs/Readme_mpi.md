# MPI for Gemini

In general Gemini uses the MPI-2 standard, which has been widely supported for over a decade by most MPI libraries.
Here's how to get MPI for common computing platforms.

## MacOS Homebrew

[Homebrew](https://brew.sh)
is a very popular development repository for MacOS.
Installing the latest MPI is simply:

```sh
brew install open-mpi
# or
brew install mpich
```

## Linux

* Ubuntu / Debian / Windows Subsystem for Linux: `apt install libopenmpi-dev openmpi-bin`
* CentOS: `dnf install openmpi-devel`

HPC users can often switch to a recent GCC version with matching MPI library.

Alternatively, [Intel oneAPI](./Linux_intel_oneapi.md)
provides Intel MPI and Scalapack on Linux.

## Windows

In general for Fortran development on Windows,
[MSYS2](https://www.scivision.dev/install-msys2-windows/)
provides a comprehensive development solution.
From the MSYS2 terminal, install MPI by:

```sh
pacman -S mingw-w64-x86_64-msmpi
```

Install
[Microsoft MS-MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi-release-notes),
which gives `mpiexec`.

Alternatively, [Intel oneAPI](./Windows_intel_oneapi.md)
provides Intel MPI and Scalapack on Windows.
We do not use MSYS2/GCC libraries with Windows oneAPI as they are ABI incompatible.
Use the oneAPI Command Prompt on Windows.

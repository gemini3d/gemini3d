# MPI for Gemini

In general Gemini uses the MPI-2 standard.
This means any MPI library from the past decade should work.
Here's how to get MPI for common computing platforms.

## MacOS Homebrew

[Homebrew](https://brew.sh)
is a very popular development repository for MacOS.
Installing the latest MPI is simply:

```sh
brew install openmpi
```

## Linux

Currently supported versions of Debian, Ubuntu, CentOS, and Windows Subsystem for Linux generally have Gfortran &ge; 6 available.

* Ubuntu / Debian / Windows Subsystem for Linux: `apt install libopenmpi-dev openmpi-bin`
* CentOS: `yum install openmpi-devel`

For CentOS, typical HPC will have the ability to switch to a recent GCC verison with matching MPI library.
If not, compile MPI--it will take about 20 minutes:

```sh
python -m gemini3d.prereqs gcc openmpi
```

The "python -m gemini3d.prereqs" command comes from
[PyGemini](https://github.com/gemini3d/pygemini)

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

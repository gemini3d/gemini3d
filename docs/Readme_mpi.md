# MPI for Gemini

In general Gemini uses the MPI-3 standard, which is widely supported by currently maintained MPI libraries.
Here's how to get MPI for common computing platforms.

## MacOS Homebrew

[Homebrew](https://brew.sh)
is a popular development repository for MacOS.
Installing the latest MPI is simply:

```sh
brew install open-mpi
# or
brew install mpich
```

Note: in general, any libraries that use MPI should be compiled with the same MPI library.

## Linux

* Ubuntu / Debian / Windows Subsystem for Linux: `apt install libopenmpi-dev openmpi-bin`
* CentOS: `dnf install openmpi-devel`

HPC users can often switch to a recent GCC version with matching MPI library.

Alternatively, [Intel oneAPI](./Linux_intel_oneapi.md)
provides Intel MPI and Scalapack on Linux.

## Windows

We suggest using Windows Subsystem for Linux from [Microsoft Store](https://apps.microsoft.com/store/detail/ubuntu-22041-lts/9PN20MSR04DW).

Alternatively,
[Intel oneAPI](./Windows_intel_oneapi.md)
provides Intel MPI and Scalapack on Windows.
Use the oneAPI Command Prompt on Windows.

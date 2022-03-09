# Build Gemini3D with GCC on Windows

If using Windows Subsystem for Linux, simply use [Linux GCC](./Linux_gcc.md) procedure.

## GCC install

[MSYS2](https://www.scivision.dev/install-msys2-windows)
provides a comprehensive Windows development solution from the Windows terminal.
From the MSYS2 terminal, install GCC C++ and Fortran compilers:

```sh
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran
```

Also install CMake, Ninja, HDF5, and MS-MPI:

```sh
pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja mingw-w64-x86_64-hdf5 mingw-w64-x86_64-msmpi
```

Install
[Microsoft MS-MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi-release-notes),
which gives `mpiexec`.

Add MSYS2 binary directory to environment variable PATH, like `C:\msys64\mingw64\bin`

## One-time Gemini3D external library setup

```sh
git clone https://github.com/gemini3d/external

cmake -S external -B external/build -DCMAKE_INSTALL_PREFIX=~/lib_gcc

cmake --build external/build
```

that installs Gemini3d external libraries under ~/lib_gcc.
This path is arbitrary but should be distinct between compilers.

## Build and Test Gemini3D

```sh
git clone https://github.com/gemini3d/gemini3d

cmake -S gemini3d -B build/gemini3d -G Ninja -DCMAKE_PREFIX_PATH=~/lib_gcc

cmake --build gemini3d/build

ctest --test-dir gemini3d/build
```

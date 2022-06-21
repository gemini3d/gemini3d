# GCC on Windows

The simplest method of using GCC is via Windows Subsystem for Linux.
[Install WSL](https://docs.microsoft.com/en-us/windows/wsl/install#install-wsl-command)
by typing in Windows Terminal:

```sh
wsl --install
```

and then use [Linux GCC](./Linux_gcc.md) procedure.

---

If you wish to use Gemini3D in native Windows instead, that is possible using MSYS2 as follows.

## MSYS2 GCC install

[MSYS2](https://www.scivision.dev/install-msys2-windows)
provides a comprehensive Windows development solution from the Windows terminal.
From the MSYS2 terminal, install CMake:

```sh
pacman -S mingw-w64-x86_64-cmake
```

Install
[Microsoft MS-MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi-release-notes),
which gives `mpiexec`.

Add MSYS2 binary directory to environment variable PATH, like `C:\msys64\mingw64\bin`

## One-time Gemini3D external library setup

```sh
git clone https://github.com/gemini3d/external

cmake -P external/scripts/requirements.cmake
# gives command to install compiler and system libraries

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

Note: to avoid having to type "-G Ninja", set environment variable `CMAKE_GENERATOR` to `Ninja`

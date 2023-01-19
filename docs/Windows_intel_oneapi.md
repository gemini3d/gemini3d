# Build Gemini3D with Intel oneAPI on Windows

Intel oneAPI (no cost) provides Intel MPI, LAPACK, and Scalapack on Windows.

Whenever wanting to use oneAPI in Windows, use the "oneAPI Command Prompt for Intel 64" on Windows.
Note: the Visual Studio CMake generator doesn't work for these project (or a lot of others).
We will use Ninja backend for CMake.

## setup Windows Intel oneAPI

Install latest no cost [Visual Studio Community](https://visualstudio.microsoft.com/vs/community/).
No particular options are needed -- a minimal install is fine.

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

## One-time Gemini3D external library setup

From the oneAPI Command prompt Window:

```sh
git clone https://github.com/gemini3d/external
```

It's necessary to use Ninja backend:

```sh
cmake -P external/scripts/install_ninja.cmake
```

Set environment variable `CMAKE_GENERATOR` to `Ninja`.
Add to environment variable `PATH` the directory containing `ninja.exe`.

Then you can build and install the external libraries:

```sh
cmake -S external -B external/build -DCMAKE_INSTALL_PREFIX=~/lib_intel

cmake --build external/build
```

that installs Gemini3d external libraries under ~/lib_intel (c:/users/username/lib_intel).
This path is arbitrary but should be distinct between compilers.

## Build and Test Gemini3D

```sh
git clone https://github.com/gemini3d/gemini3d

cmake -S gemini3d -B build/gemini3d -G Ninja -DCMAKE_PREFIX_PATH=~/lib_intel

cmake --build gemini3d/build

ctest --test-dir gemini3d/build
```

Note: to avoid having to type "-G Ninja", set environment variable `CMAKE_GENERATOR` to `Ninja`

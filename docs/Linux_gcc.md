# Build Gemini3D with GCC on Linux

This method also works for
[Windows WSL](https://docs.microsoft.com/en-us/windows/wsl/install).

GCC 7.5 and newer works with Gemini3D.

```sh
apt install cmake
# or
dnf install cmake
```

## One-time Gemini3D external library setup

```sh
git clone https://github.com/gemini3d/external

cmake -P external/scripts/requirements.cmake
# gives command to install compiler and system libraries

cmake -P external/build-online.cmake
```

that installs Gemini3d external libraries under ~/libgem_gcc.
This path is arbitrary but should be distinct between compilers.

NOTE: If CMake is too old, install a new CMake:

```sh
cmake -P external/scripts/install_cmake.cmake
```

## Build and Test Gemini3D

```sh
git clone https://github.com/gemini3d/gemini3d

cmake -S gemini3d -B build/gemini3d -DCMAKE_PREFIX_PATH=~/libgem_gcc

cmake --build gemini3d/build

ctest --test-dir gemini3d/build
```

## Troubleshooting

The compiler relies on libc and libstdc++.
If build errors about "filesystem" at C++ link time, the system configuration may be messed up.
If problems, ensure environment variable LD_LIBRARY_PATH has first the libc and libstdc++ for the GCC version you wish to use.

# Build Gemini3D with GCC on Linux

GCC 7.5 and newer works with Gemini3D.

```sh
apt install g++ gfortran
# or
dnf install gcc-c++ gcc-gfortran
```

The compiler relies on libc and libstdc++.
If build errors about "filesystem" at C++ link time, the system configuration may be messed up.
If problems, ensure environment variable LD_LIBRARY_PATH has first the libc and libstdc++ for the GCC version you wish to use.

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

# Build Gemini3D with GCC on MacOS

AppleClang LLVM compilers work fine with GCC/Gfortran in general on MacOS.

Install the latest GCC/Gfortran via [Homebrew](https://brew.sh):

```sh
brew install gcc
```

The compiler relies on libc and libstdc++.
If build errors missing "stdio.h" or "filesystem" at C++ link time, the system configuration may be messed up.

If problems, try installing Xcode:

```sh
xcode-select --install
```

and try these environment variables:

```sh
export LIBRARY_PATH=$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
export CPLUS_INCLUDE_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include

export CXXFLAGS=-I$CPLUS_INCLUDE_PATH
export CFLAGS=$CXXFLAGS
```

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

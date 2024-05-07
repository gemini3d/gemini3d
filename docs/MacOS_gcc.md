# Build Gemini3D with GCC on MacOS

GCC compilers are recommended for macOS as AppleClang conflicts with Gfortran for exception handling.

Install libraries via [Homebrew](https://brew.sh) such as CMake:

```sh
brew install cmake
```
## Build and Test Gemini3D

```sh
git clone https://github.com/gemini3d/gemini3d

cmake -S gemini3d -B build/gemini3d

cmake --build gemini3d/build

ctest --test-dir gemini3d/build
```

## Troubleshooting

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

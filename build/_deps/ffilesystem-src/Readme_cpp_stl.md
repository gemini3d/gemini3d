# C++ standard library compatibility

Ffilesystem relies on the C++ standard library (STL) for its implementation as do most C++ programs.
STL is a critical part of the C++ programming environment.
Computing platforms provide STL in distinct ways:

* Windows MSVC: STL is provided by the [Visual Studio C++ compiler](https://github.com/microsoft/STL). [Changelog](https://github.com/microsoft/STL/wiki/Changelog)
* Windows [MSYS2](https://www.msys2.org/docs/environments): STL is provided by [GNU libstdc++](https://packages.msys2.org/packages/mingw-w64-ucrt-x86_64-gcc) or [LLVM libc++](https://packages.msys2.org/packages/mingw-w64-ucrt-x86_64-libc++).
* macOS AppleClang: STL is provided by [Xcode](https://developer.apple.com/forums/thread/715385) and is not a specific file.
* macOS [Homebrew](https://brew.sh): STL is by default the macOS Apple system STL, but can be specified to be Homebrew [LLVM libc++](https://formulae.brew.sh/formula/llvm) or [GNU libstdc++](https://formulae.brew.sh/formula/gcc).

Where allowed by the complier/linker/libc++, one can specify the libc++ library like:

```sh
cmake -B build -DCMAKE_CXX_FLAGS="-stdlib=libc++"
```


On Linux systems, STL is provided by default by the GNU libstdc++ library for many (almost all popular) distributions.
Other Linux distributions such as
[Chimera Linux](https://chimera-linux.org/)
use
[LLVM libc++](https://archive.fosdem.org/2024/schedule/event/fosdem-2024-2555-building-a-linux-distro-with-llvm).

The STL version might need additional flags to be specified in the build system.
For example, on Linux with Intel oneAPI or NVHPC compilers, libstdc++ is the default STL.
If the STL version is older, linker flags such as `-lstdc++fs -lstdc++` might be needed.
A good way to detect this need is by checking the STL version by CMake `try_run` or similar, and applying the flags conditionally.
This is done in Ffilesystem cmake/compilers.cmake script.

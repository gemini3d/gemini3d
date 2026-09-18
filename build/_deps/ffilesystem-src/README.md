# Ffilesystem: platform-independent, compiler-agnostic path manipulation and filesystem library

[![DOI](https://zenodo.org/badge/433875623.svg)](https://zenodo.org/badge/latestdoi/433875623)
[![ci](https://github.com/scivision/ffilesystem/actions/workflows/ci.yml/badge.svg)](https://github.com/scivision/ffilesystem/actions/workflows/ci.yml)
[![ci_windows](https://github.com/scivision/ffilesystem/actions/workflows/ci_windows.yml/badge.svg)](https://github.com/scivision/ffilesystem/actions/workflows/ci_windows.yml)
[![oneapi-linux](https://github.com/scivision/ffilesystem/actions/workflows/oneapi-linux.yml/badge.svg)](https://github.com/scivision/ffilesystem/actions/workflows/oneapi-linux.yml)
[![ci_fpm](https://github.com/scivision/ffilesystem/actions/workflows/ci_fpm.yml/badge.svg)](https://github.com/scivision/ffilesystem/actions/workflows/ci_fpm.yml)

Platform independent (Linux, macOS, Windows, Android, Cygwin, WSL, BSD, ...) and compiler-agnostic Ffilesystem path manipulation library.
Simplicity and efficiency are focuses of Ffilesystem.
Ffilesystem backend is implemented in C++17, with C++ functions generally using
[<string_view>](https://en.cppreference.com/w/cpp/header/string_view.html)
inputs for simplicity and speed.
C and Fortran interfaces are provided for nearly all C++ backend functions.
For the Fortran interface, native Fortran `character` type is used.
For the C interface,
[C-style null-terminated strings](https://en.cppreference.com/w/cpp/string/byte.html)
are used.
If available,
[C++ standard library `<filesystem>`](https://en.cppreference.com/w/cpp/filesystem)
is used.
The C++ backend accesses the
[C standard library](https://en.wikipedia.org/wiki/C_standard_library)
to access filesystem and system parameters.
Networked file systems and FUSE (e.g. SSHFS) are supported as well as local filesystems.

Ffilesystem supports
[UTF-8](https://utf8everywhere.org/)
via C++ `std::string` or Fortran `character` as the path input / output argument types.

Ffilesystem does not throw or catch C++ exceptions itself.

Ffilesystem header
[ffilesystem.h](./include/ffilesystem.h)
can be used from C and C++ project code--see
[example](./example).
The C interface allows reuse of Ffilesystem functions in other code languages such as
[Matlab](./matlab/Readme.md)
via FFI.

The optional Fortran interface is built by default.
Disable Fortran by

```sh
cmake -Dffilesystem_fortran=false -Bbuild
```

Ffilesystem brings full, fast filesystem functionality to Fortran.

The language standards must be at least:

* C++17 standard library [STL](./Readme_cpp_stl.md)
* (optional) Fortran 2003

Ffilesystem works with popular C++ STL and C standard library implementations including:
[glibc](https://sourceware.org/glibc/),
[newlib](https://sourceware.org/newlib/),
[musl](https://musl.libc.org/),
[Android LLVM libc++](./Readme_android.md),
[iPhone / iOS builds](./Readme_ios.md),
[Cosmopolitan universal binaries](./Readme_cosmopolitan.md),
[macOS universal binaries](./Readme_macos.md),
BSD libc,
[Microsoft CRT](https://en.wikipedia.org/wiki/Microsoft_Windows_library_files#CRT),
among others.
On Linux and Cygwin, symbol
[_DEFAULT_SOURCE](https://man7.org/linux/man-pages/man7/feature_test_macros.7.html)
is defined if needed to enable POSIX functions that may not be defined in the C standard library if ANSI `-std=c++*` is used instead of `-std=gnu++*`.
RAII std::string buffers are used for all string representations of paths including for C API and system calls, and are automatically freed.
Ffilesystem
[std::string buffers](https://learn.microsoft.com/en-us/archive/msdn-magazine/2015/july/c-using-stl-strings-at-win32-api-boundaries)
are dynamically sized according to the actual path length.
See [Readme_memory](./Readme_memory.md)

For Windows drive letters without a slash after the colon, the path is treated as
[a relative path](https://learn.microsoft.com/en-us/windows/win32/fileio/naming-a-file).
This is how the Windows API, ComSpec, and C++ STL work.

Inspired by other languages' filesystem libraries, Ffilesystem provides a comprehensive set of functions for path manipulation and filesystem operations.
The API is designed to be familiar to users of other languages' filesystem libraries, such as:

* [Python pathlib](https://docs.python.org/3/library/pathlib.html)
* [Julia filesystem](https://docs.julialang.org/en/v1/base/file/)
* [Go io/fs](https://pkg.go.dev/io/fs)
* [Rust std::fs](https://doc.rust-lang.org/std/fs/index.html)
* [C# / .NET System.IO](https://learn.microsoft.com/en-us/dotnet/api/system.io)
* [Ruby File](https://ruby-doc.org/3.4.1/File.html)
* [Node.js fs](https://nodejs.org/api/fs.html)
* [Java NIO Files](https://docs.oracle.com/en/java/javase/25/docs/api/java.base/java/nio/file/package-summary.html)
* [stdlib for Matlab](https://github.com/geospace-code/matlab-stdlib)
* [C++ STL filesystem](https://en.cppreference.com/w/cpp/filesystem)
* [C standard library](https://en.cppreference.com/w/c/io)
Important Ffilesystem functions are
[benchmarked](./benchmark.md)
to help improve performance.
Advanced / conceptual development takes place in
[ffilesystem-concepts](https://github.com/scivision/ffilesystem-concepts).

## Compiler support

Ffilesystem supports compilers including:

* GCC &ge; 7 (gcc/g++, gfortran)
* LLVM Clang &ge; 9 (clang/clang++, flang or gfortran)
* Intel oneAPI &ge; 2023.1 (icx, icpx, ifx)
* Intel Classic &ge; 2021.9 (icpc, ifort)
* AMD AOCC (clang/clang++, flang)
* NVIDIA HPC SDK (nvc++, nvfortran)
* Visual Studio (C/C++)
* Cray: using Cray compilers alone (cc, CC, ftn) or using GCC or Intel backend

To help debug with older compilers, disable C++ `<filesystem>`:

```sh
cmake -Bbuild -Dffilesystem_cpp=off
```

### libstdc++

Some systems have broken, obsolete, or incompatible libstdc++.

**Intel**: If Intel compiler linker errors use GCC >= 9.1.
This can be done by setting environment variable CXXFLAGS to the top level GCC >= 9.1 directory.
Set environment variable CXXFLAGS for
[Intel GCC toolchain](https://www.intel.com/content/www/us/en/develop/documentation/oneapi-dpcpp-cpp-compiler-dev-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/compatibility-options/gcc-toolchain.html)
like:

```sh
export CXXFLAGS=--gcc-toolchain=/opt/rh/gcc-toolset-10/root/usr/
```

which can be determined like:

```sh
scl enable gcc-toolset-10 "which g++"
```

## Build

Ffilesystem can be built with CMake or Fortran Package Manager (FPM).

"libffilesystem.a" is the library binary built that contains the Fortran "filesystem" module--it is the only binary you need to use in your project.

Please see the [API docs](./API.md) for extensive list of functions/subroutines.

Use any one of these methods to build Ffilesystem.
The self-tests are optional and not built by default.
The tests use Boost.UT framework.

### CMake

```sh
cmake -B build
cmake --build build
```


Optionally, build and run the self-tests:

```sh
cmake -B build -Dffilesystem_BUILD_TESTING=on
cmake --build build
ctest --test-dir build
```

The default library with CMake is static; to build shared library:

```sh
cmake -B build -DBUILD_SHARED_LIBS=on
...
```

### Fortran Package Manager (FPM):

```sh
fpm --cxx-flag=-std=c++17 build
# c++17 is the minimum, can use newer
```

---

GNU Make:

```sh
make
```

We provide Fortran REPL "filesystem_cli" and C++ REPL "fs_cli" for interactive testing of Ffilesystem routines.

### Build options

We have begun to implement Fortran standard
[ISO_Fortran_binding.h (see section 18.5)](https://j3-fortran.org/doc/year/18/18-007r1.pdf)
Fortran bindings for C++ functions.
This is by default OFF and can be enabled with CMake option `cmake -Dffilesystem_iso_fortran_binding=on`.
A CMake workflow that tests the ISO_Fortran_binding.h Fortran bindings is:

```sh
cmake --workflow isofortran
```

---

Fortran "filesystem" module contains OPTIONAL (enabled by default) Fortran type "path_t" that contains properties and methods.
The "path_t" type uses getter and setter procedure to access the path as a string `character(:), allocatable`.

```fortran
use filesystem, only : path_t

type(path_t) :: p

p = path_t("my/path")  !< setter

print *, "path: ", p%path() !< getter
```

CMake detects if Fortran 2003 `type` is available and enables `path_t` by default.
To manually enable / disable `path_t` with CMake set command option `cmake -DHAVE_F03TYPE=1` or `cmake -DHAVE_F03TYPE=0` respectively.

---

[statx()](https://www.man7.org/linux/man-pages/man2/statx.2.html)
is used by default on
[glibc](https://www.gnu.org/software/gnulib/manual/html_node/statx.html)
version
[2.28 or newer](https://www.sourceware.org/glibc/wiki/Release/2.28)
systems to get file information.
Like
[Python](https://github.com/python/cpython/pulls?q=is%3Apr+statx+)
[os.statx()](https://docs.python.org/3.15/library/os.html#os.statx),
we check for `statx()` availability at compile time on Linux and Android and fall back at runtime to `stat()`.

## Self test

The optional self-tests provide reasonable coverage of the Ffilesystem library.
Several of the tests use `argv[0]` as a test file.
We are aware of the shortcomings of `argv[0]` to get the executable name.
We provide the function `fs_exepath()` to get the executable path reliably.

## Usage from other projects

The [example](./example) directory contains a use pattern from external projects.
One can either `cmake --install build` ffilesystem or use CMake ExternalProject or
[FetchContent](https://gist.github.com/scivision/245a1f32498d15a87a507051857327d9)
from the other project.
To find ffilesystem in your CMake project:

```cmake
find_package(ffilesystem CONFIG REQUIRED)
```

CMake package variables `ffilesystem_cpp` and `ffilesystem_fortran` indicate whether ffilesystem was built with C++ `<filesystem>` and/or Fortran support.

[ffilesystem.cmake](./cmake/ffilesystem.cmake)
would be included from the other project to find or build Ffilesystem automatically.
It provides the appropriate imported targets for shared or static builds, including Windows DLL handling.

## Notes

GCC 6.x and older aren't supported due to lack of C++17 string_view support.

In C++ code, we generally use the
[size_type](https://en.cppreference.com/w/cpp/types/size_t)
of the class as a best practice--for example `std::string::size_type` where appropriate instead of `std::size_t`.
We use `size_t` at the C interfaces for clarity and also certain internal library calls.
`ssize_t` is used in certain non-Windows internal-only function calls.

For those desiring to use Windows compiler ClangCL that comes with Visual Studio, this is enabled by CMake preset or workflow "clangcl".
Fortran FlangCL compiler isn't yet working as iso_c_binding is missing and CMake may need support added.
The computer needs Visual Studio Installer "Individual Components":

* "C++ Clang Compiler for Windows"
* "MSBuild support for LLVM (clang-cl) toolset"

### Possible future features

* inquire if a file is encrypted or compressed, etc.
* determine if the terminal is in a [VT100 compatible mode](https://gitlab.kitware.com/cmake/cmake/-/merge_requests/10759/diffs)

### Other C++ filesystem libraries

Ffilesystem emphasizes simplicity and reasonable performance and reliability for scientific computing, particularly on HPC systems.
A highly performance-oriented C++ low-level no TOCTOU filesystem library is
[LLFIO](https://github.com/ned14/llfio).
An older C++ object-oriented interface is
[CppFS](https://github.com/cginternals/cppfs).

Other implementations of C++ filesystem include:

* [Boost.Filesystem](https://www.boost.org/doc/libs/1_86_0/libs/filesystem/doc/index.htm) what the stdlib filesystem is based on, often tries newer features. [Boost.Filesystem source code](https://github.com/boostorg/filesystem)
* [ghc-filesystem](https://github.com/gulrak/filesystem) for older compilers.
* deprecated `<experimental/filesystem>` is missing vital lexical operations.

### Other Fortran filesystem libraries

Other Fortran libraries that provide interfaces to filesystems include the following.
Generally they have noticeably fewer functions than Ffilesystem.
They typically implement many functions in Fortran, where with Ffilesystem we implement in C++ or C++ `<filesystem>` if available.
Ffilesystem Fortran code is optional, and is just a thin wrapper around the C++ functions.

* [stdlib_os](https://github.com/MarDiehl/stdlib_os)
* [fortyxima](https://bitbucket.org/aradi/fortyxima/src/develop/)
* [Fortran-stdlib](https://github.com/fortran-lang/stdlib/issues/201)
* [M_system](https://github.com/urbanjost/M_system) focuses on interfaces to libc

### features we won't add

There is no "is_musl()" function due to MUSL designers
[not providing](https://stackoverflow.com/questions/58177815/how-to-actually-detect-musl-libc)
a
[MUSL feature macro](https://wiki.musl-libc.org/faq.html).

Disk / partition formatting is something we won't add. We do have functions to report the partition type and available capacity of the disk partition.

### Windows

Like
[Microsoft STL](https://github.com/microsoft/STL/issues/2256),
Ffilesystem is not designed for UNC "long" paths.
We recommend using a UNC path to a mapped drive letter.
Windows
[long paths](https://github.com/microsoft/STL/issues/1921)
are partially implemented due to
[limitations](https://www.boost.org/doc/libs/1_88_0/libs/filesystem/doc/reference.html#windows-path-prefixes).
Boost.Filesystem handles long paths and we may implement Boost.Filesystem as an optional backend in the future.

Enable Windows
[developer mode](https://learn.microsoft.com/en-us/windows/apps/get-started/developer-mode-features-and-debugging)
to use symbolic links if needed.

### C++ filesystem discussion

Security
[research](https://www.reddit.com/r/cpp/comments/151cnlc/a_safety_culture_and_c_we_need_to_talk_about/?rdt=62365)
led to
[TOCTOU](https://en.wikipedia.org/wiki/Time-of-check_to_time-of-use)-related
patches to the C++ filesystem library in various C++ standard library implementations noted in that discussion.
Ffilesystem does NOT use remove_all, which was the TOCTOU concern addressed above.

Since the underlying C++ filesystem is not thread-safe, race conditions can occur if multiple threads are accessing the same filesystem object regardless of the code language used in the Ffilesystem library.

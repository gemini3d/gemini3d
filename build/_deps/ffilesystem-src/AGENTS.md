# Ffilesystem fast filesystem library

## Project overview

This is a lightweight C++17 or newer library for filesystem path manipulation and querying.
It is designed to work with NUL-terminated strings--for brevity and speed we always assume any string input is valid and NUL-terminated.
We assume the compiler is at least C++17 capable.
We use std::string_view heavily for speed and to avoid unnecessary allocations and copies.
We do a CMake (or Meson) configure-time check to ensure C++17 `<filesystem>` is available and linkable on the platform, falling back to our own implementations if not. The macro `HAVE_CXX_FILESYSTEM` is defined if `<filesystem>` is available and working.

## Platform support

Target platforms are Linux, Windows, Windows Subsystem for Linux (WSL), macOS, Android, and BSD.

* We aim to support a range of Linux C and C++ standard libraries--for example, the libc and libstdc++ available with currently supported Ubuntu LTS and RHEL versions.
* We aim to support relative recent BSD flavors (FreeBSD, OpenBSD, NetBSD).
* We aim to support versions of macOS currently supported by Apple.
* We aim to support all freely available, currently maintained compiler versions, including GCC/GFortran, Clang/Flang, Intel oneAPI, and NVIDIA HPC SDK.

### Windows support

We support Windows 10 and newer -- whatever is currently supported by Microsoft.

Compiler: support MinGW (GCC, Clang), MSVC, and Intel oneAPI.

UNC paths are generally not supported--with MSVC and GCC, `std::filesystem` itself does not support UNC paths.
Only Boost::Filesystem supports UNC paths, but we haven't yet implemented Boost::Filesystem .

## Function API

Generally, the paths are passed in as `std::string_view` and returned as `std::string` or `std::optional<std::string>` if the operation can fail.
We return `std::string` instead of `std::string_view` because the C API is important and generally the users are not interested in managing the lifetime of the string_view.

* fs_print_error() is polymorphic and grabs GetLastError() on Windows or errno on POSIX systems. Assume all arguments to fs_print_error() are of the correct type.

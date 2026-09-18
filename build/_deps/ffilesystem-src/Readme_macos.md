# macOS Universal binaries

[macOS universal binaries](https://en.wikipedia.org/wiki/Universal_binary)
are supported by Ffilesystem.
To build a universal binary with CMake:

```sh
cmake -DCMAKE_OSX_ARCHITECTURES="arm64;x86_64" -B build -Dffilesystem_fortran=off

cmake --build build -t fs_cli
```

To force the universal binary to run in x86_64 mode from an arm64 Mac and verify Rosetta 2 is being used:

```sh
arch -x86_64 build/app/fs_cli

Ffs> is_rosetta
1

Ffs> cpu_arch
x86_64
```

Whereas using the default native mode on an arm64 Mac:

```sh
build/app/fs_cli

Ffs> is_rosetta
0

Ffs> cpu_arch
arm64
```

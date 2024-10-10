# CMake options

Note: a common convention is to compile (build) project in the "build/" directory.
That is, the Gemini.bin executable and other libraries and test executables will be created under gemini3d/build.
Feel free to use any directory name you like.

---

When using a meta-build system like CMake, building a program takes two steps.

1. configure: "cmake -B build": user specifies options (if any) and CMake determines what compiler options and finds libraries.
2. build: "cmake --build build": CMake uses a low-level build system, typically Ninja or Make, to orchestrate the C, C++, and Fortran compiler commands in parallel.

If a user makes changes to source code, they need only rerun the `cmake --build build` command.
If a user wishes to change a CMake option, run both commands again.

## Gemini3D build

Many Gemini3D users don't need any custom options to get started.
In that case, build Gemini3D from the gemini3d/ directory by:

```sh
cmake -B build

cmake --build build
```

## CMake configure

CMake build-time options in general are enabled or disabled at *configure time* like

```sh
cmake -B build -Doption=true

cmake -B build -Doption=false
```

If you've already built Gemini but wish to change a CMake configuration option, you can reconfigure as above, then rebuild:

```sh
cmake --build build
```

Those adding or modifying Gemini3D code itself may be interested in
[Debug builds](./Readme_debug.md).

### macOS

macOS can use the default AppleClang compilers with Gfortran.
The GCC / Gfortran compilers are available from Homebrew, Macports, etc.

When specifying GCC, macOS itself provides a "fake" `gcc` that is linked to Clang.
In general, to specify the actual GNU GCC, set environment variables like:

```sh
export FC=gfortran-14 CC=gcc-14 CXX=g++-14
```

### GLOW

NCAR GLOW is optional.
Auroral emissions use GLOW.

Disable GLOW by:

```sh
cmake -B build -Dglow=off
```

### MSIS 2.x

The neutral atmosphere model MSISE00 is used by default.
To use the newer MSIS 2.x, the simulation config.nml must specify the following to actually use MSIS 2.x:

```ini
&neutral_BG
msis_version = 21  ! MSIS version multiplied by 10 e.g. 21 is MSIS 2.1
/
```

Omitting this namelist variable or specifying `msis_version=0` uses MSISE00.

### HWM14

Gemini3D may use the HWM14 horizontal wind model by:

```sh
cmake -B build -Dhwm14=on
```

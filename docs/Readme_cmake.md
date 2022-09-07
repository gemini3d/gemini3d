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

## One-time Gemini3D external library setup

```sh
git clone https://github.com/gemini3d/external

cmake -P external/scripts/requirements.cmake
# gives command to install compiler and system libraries

cmake -S external -B external/build -DCMAKE_INSTALL_PREFIX=~/libgem
# ~/libgem is an arbitrary location--pick whatever name you like.

cmake --build external/build
```

that installs Gemini3d external libraries under ~/libgem.

## Gemini3D build

Many Gemini3D users don't need any custom options to get started.
In that case, build Gemini3D from the gemini3d/ directory by:

```sh
cmake -B build -DCMAKE_PREFIX_PATH=~/libgem
# ~/libgem is the arbitrary location you installed Gemini3D/external libraries to.

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

### MacOS

MacOS users usually use:

* Homebrew GCC
* Macports GCC
* Intel oneAPI

When trying to use GCC, be aware MacOS itself provides a "fake" `gcc` that is linked to Clang.
You may desire to use GCC instead of Clang if you get compiler errors.
These errors may arise due to certain ABI incompatibilities between Clang and GCC.

Tell CMake (or other build systems) you want to use just GCC by:

```sh
export FC=gfortran-11 CC=gcc-11 CXX=g++-11
```

To default to GCC / Gfortran, add the line above in ~/.zprofile.

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

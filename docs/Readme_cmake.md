# CMake options

Note: a common convention is to compile (build) project in the "build/" directory.
That is, the Gemini.bin executable and other libraries and test executables will be created under gemini3d/build.
Feel free to use a different directory name if you like.

---

When using a meta-build system like CMake, building a program takes two steps.

1. configure: user specifies options (if any) and CMake determines what compiler options and finds libraries. If libraries are missing, CMake prepares to build them.
2. build: CMake uses a low-level build system, typically Ninja or GNU Make, to orchestrate the C and Fortran compiler commands in parallel.

If a user makes changes to source code, they need only rerun the build command.
If a user wishes to change a CMake option, run both commands again.

## Normal Gemini build

Many Gemini3D users do not need any custom options.
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

### Windows

If you happen to look inside build/CMakeFiles/CMakeError.log, with MS-MPI you will see text like

```
build/CMakeFiles/FindMPI/test_mpi.f90:2:11:

    2 |       use mpi_f08
      |           1
Fatal Error: Cannot open module file 'mpi_f08.mod' for reading at (1): No such file or directory
compilation terminated.
```

This is normal because MS-MPI does not yet have MPI-3.
Gemini3D uses MPI-2 so this is not relevant to Gemini3D.

### MacOS: Homebrew

If using Homebrew on MacOS, be sure Homebrew's GCC is used instead of AppleClang or other non-Homebrew compilers so that the Homebrew library ABIs match the compiler ABI.

```sh
FC=gfortran-10 CC=gcc-10 cmake -B build
```

If you always use GCC / Gfortran, set the environment variables in ~/.bashrc or ~/.zshenv like:

```sh
export FC=gfortran-10 CC=gcc-10
```

## GLOW

NCAR GLOW is automatically installed, but optional in general.
Auroral emissions use GLOW.

Disable GLOW by:

```sh
cmake -B build -Dglow=off
```

## HDF5

HDF5 is enabled by default, and disabled by:

```sh
cmake -B build -Dhdf5=off
```

If desired, build HDF5 manually by
[PyGemini](https://github.com/gemini3d/pygemini):

```sh
python -m gemini3d.prereqs gcc hdf5
```

## NetCDF

NetCDF is disabled by default, and enabled by:

```sh
cmake -B build -Dnetcdf=on
```

## MSIS 2.0

The neutral atmosphere model MSISE00 is used by default.
To use the newer MSIS 2.0:

```sh
cmake -B build -Dmsis20=on
```

Additionally, the simulation config.nml must specify the following to actually use MSIS 2.0:

```ini
&neutral_BG
msis_version = 20
/
```

Omitting this namelist variable or specifying `msis_version=0` uses MSISE00.

## HWM14

Gemini3D may use the HWM14 horizontal wind model by:

```sh
cmake -B build -Dhwm14=on
```

# build-time options

CMake build-time options in general are enabled or disabled like

```sh
cmake -B build -Doption=true

cmake -B build -Doption=false
```

## Homebrew

If using Homebrew on MacOS, be sure Homebrew's GCC is used instead of AppleClang or other non-Homebrew compilers so that the Homebrew library ABIs match the compiler ABI.

```sh
FC=gfortran-9 CC=gcc-9 cmake -B build

cmake --build build
```

If you always use GCC / Gfortran, set the environment variables in ~/.bashrc or ~/.zshenv like:

```sh
export FC=gfortran-9
export CC=gcc-9
```

## GLOW

NCAR GLOW is automatically installed, but optional in general.
Auroral emissions use GLOW.

Disable GLOW by:

```sh
cmake -B build -Dglow=off

cmake --build build
```

## HDF5

HDF5 is enabled by default, and disabled by:

```sh
cmake -B build -Dhdf5=off

cmake --build build
```

If you would like to build HDF5 yourself instead of installing it via your package manager, you can use CMake from the top gemini3d/ directory:

```sh
ctest -S scripts/build_hdf5.cmake
```

or by PyGemini (assuming it was setup during your Gemini3D setup or manually):

```sh
gemini_prereqs gcc hdf5
```

## NetCDF

NetCDF is disabled by default, and enabled by:

```sh
cmake -B build -Dnetcdf=on

cmake --build build
```

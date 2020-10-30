# build-time options

CMake build-time options in general are enabled or disabled like

```sh
cmake -B build -Doption=true

cmake -B build -Doption=false
```

### Windows

If you happen to look inside build/CMakeFiles/CMakeError.log, with MS-MPI you will currently see text like

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

If you would like to build HDF5 yourself instead of installing it via your package manager, by PyGemini (assuming it was setup during your Gemini3D setup or manually):

```sh
python -m gemini3d.prereqs gcc hdf5
```

## NetCDF

NetCDF is disabled by default, and enabled by:

```sh
cmake -B build -Dnetcdf=on

cmake --build build
```

# Installing GEMINI Prerequisites

As indicated in the main readme there are a number of prerequisites for GEMINI.  This readme covers some details of how to install these on various platforms.  We focus here on ease of use and describe the "least effort" methods for getting needed software.

## Working with MacOS

You will need to install XCode through the app store.
Then it may be necessary (depending on your OS version) to manually install XCode command line tools:

```sh
xcode-select --install
```

It is strongly recommended that you install Homebrew by following instructions on [the homebrew website](https://brew.sh).
Instructions in this readme assume that you will use this package manager to get most of the prereqs.  Users have also reported that [Macports](https://www.macports.org) works fine for getting required packages, as well.  Both are available from the linked websites.

## Installing Compilers

Many default installations will not have the required compilers, e.g. there is no default Fortran compiler in Mac OS, and many Linux distributions install without a C++ compiler.
See [Readme_compilers](./Readme_compilers.md) for more info.

## Installing Parallelization Libraries

MPI is often installed using a package manager since it can take a very long time to build from source.
See [Readme_mpi](./Readme_mpi.md) for more info.

## Python

Python is required for a number of GEMINI operations.
We recommend installing [miniconda](https://docs.conda.io/en/latest/miniconda.html).
Open source distributions are available for download at these websites.
As a last resort, [build Python from scratch](https://github.com/gemini3d/cmake-python-build).

## MATLAB

Extensive scripting front-ends for simulation preparation and analysis exist in MATLAB.
These are not required (implementations exist in python); however they are very useful and in some cases contain functionality not yet implemented in python.
MATLAB can be obtained from the [Mathworks web site](https://www.mathworks.com) and requires a paid license for use.

## CMake

[CMake](./Readme_cmake_install.md) is required.

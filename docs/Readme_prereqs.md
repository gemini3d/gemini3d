# Installing GEMINI Prerequisites

As indicated in the main readme there are a number of prerequisites for GEMINI.  This readme covers some details of how to install these on various platforms.  We focus here on ease of use and describe the "least effort" methods for getting needed software.  

## Working with MacOS

You will need to install XCode through the app store.  Then it may be necessary (depending on your OS version) to manually install XCode command line tools:

```sh
xcode-select --install
```

It is strongly recommended that you install homebrew by following instructions on [the homebrew website](https://brew.sh).  Instructions in this readme assume that you will use this package manager to get most of the prereqs.  Users have also reported that [Macports](https://www.macports.org) works fine for getting required packages, as well.  Both are available from the linked websites.  

# Installing Compilers

Many default installations will not have the required compilers, e.g. there is no default fortran compiler in Mac OS, and many ubuntu distributions install without a C++ compiler.  

## MacOS

To install the gcc compiler suite:

```sh
brew install gcc
```

This will also give you gfortran, by default C/C++ will be compiled with AppleClang which seems to work fine for all the installations we have worked with.  

## Linux

```sh
apt-get install gfortran
apt-get install g++
```


# Installing Parallelization Libraries

We most often use openmpi which can/should be installed using a package manager since it can take a very long time to build from source.  

## MacOS

```sh
brew install openmpi
```

## Linux

```sh
apt-get install openmpi openmpi-dev
```

# Python 3 Installation

Python 3 is required for a number of GEMINI operations; we recommend installing either [anaconda](https://www.anaconda.com) or [miniconda](https://docs.conda.io/en/latest/miniconda.html).  Open source distributions are available for download at these websites.  

# MATLAB Installation

Extensive scripting front-ends for simulation preparation and analysis exist in MATLAB.  These are not required (implementations exist in python); however they are very useful and in some cases contain functionality not yet implemented in python.  MATLAB can be obtained from the [Mathworks web site](https://www.mathworks.com) and requires a paid license for use.  

# Installing cmake

Cmake is used extensively for GEMINI's build system and its continuous integration system and so is required for any deployments.  

## MacOS

On MacOS a reasonably updated version of cmake can be obtained from your package manager, e.g.:

```sh
brew install cmake
```

## Linux

Most linux distributions (even recent ones) come with an outdated cmake package, and many do not install it by default.  There are essentially two easy ways to get an updated cmake:

1. Install the distribution cmake package (which will be outdated) and use it with our update script to get a newer version, e.g.:

```sh
apt-get install cmake
cd <GEMINI source code directory>
cmake -P scripts/install_cmake.cmake
```

2. Install cmake from source by first downloading and extracting the latest version from [cmake.org](https://cmake.org):

```sh
cd <cmake source code directory>
./bootstrap --prefix=<install directory>
make -j
make install
```
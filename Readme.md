# GEMINI

[![DOI](https://zenodo.org/badge/146920930.svg)](https://zenodo.org/badge/latestdoi/146920930)
[![CDash](./tests/data/CDash.png)](https://my.cdash.org/index.php?project=Gemini3D)
![ci_build](https://github.com/gemini3d/gemini/workflows/ci_build/badge.svg)
![ci_linux](https://github.com/gemini3d/gemini/workflows/ci_linux/badge.svg)
![ci_macos](https://github.com/gemini3d/gemini/workflows/ci_macos/badge.svg)
![ci_windows](https://github.com/gemini3d/gemini/workflows/ci_windows/badge.svg)

The GEMINI model (*G*eospace *E*nvironment *M*odel of *I*on-*N*eutral *I*nteractions) is a three-dimensional ionospheric fluid-electrodynamic model.
GEMINI is used for various scientific studies including:

* effects of auroras on the terrestrial ionosphere
* natural hazard effects on the space environment
* effects of ionospheric fluid instabilities on radio propagation

The detailed mathematical formulation of GEMINI is included in
[GEMINI-docs](https://github.com/gemini3d/GEMINI-docs).
A subroutine-level set of
[inline generated](./docs/Readme_docs.md)
documentation describing functions of individual program units is given via source code comments which are
[rendered as webpages](https://gemini3d.github.io/GEMINI/).
GEMINI uses generalized orthogonal curvilinear coordinates and has been tested with dipole and Cartesian coordinates.
Please open a
[GitHub Issue](https://github.com/gemini3d/gemini/issues)
if you experience difficulty with GEMINI.

Generally, the Git `master` branch has the current development version and is the best place to start, while more thoroughly-tested releases happen regularly.
Specific releases corresponding to published results are generally noted in the corresponding journal article.

## Platform agnostic

Gemini is OS / CPU arch / platform / compiler agnostic.
Operating system support includes Linux, MacOS, and Windows.
CPU arch support includes Intel / AMD, ARM and IBM POWER.
GEMINI can run on hardware ranging from a Raspberry Pi to laptop to a high-performance computing (HPC) cluster.
One could run large 2D or very small 3D simulations (not exceeding a few million grid points) on a quad-core workstation, but may take quite a while to complete.

For large 3D simulations (more than 20M grid points), GEMINI should be run in a cluster environment or a "large" multi-core workstation (e.g. 12 or more cores).
Runtime depends heavily on the grid spacing used, which determines the time step needed to insure stability,
For example we have found that a 20M grid point simulations takes about  4 hours on 72 Xeon E5 cores.  200M grid point simulations can take up to a week on 256 cores.
It has generally been found that acceptable performance requires > 1GB memory per core; moreover, a large amount of storage (hundreds of GB to several TB) is needed to store results from large simulations.

## Quick start

To build Gemini and run self-tests takes about 10 minutes on a laptop.

Requirements:

* Fortran 2008 compiler. See [compiler help](./docs/Readme_compilers.md) if needed.
  * Gfortran / GCC &ge; 7
  * Intel oneAPI HPC Toolkit (free to use for all)
* [CMake](https://cmake.org/download/): if your CMake is too old, update by running `cmake -P scripts/install_cmake.cmake` or from Python `pip install cmake`
* Git: the Gemini3D software stack uses Git to version lock reproducible builds.

Recommended:

* MPI: any of OpenMPI, IntelMPI, MPICH, MS-MPI. See [MPI help](./docs/Readme_mpi.md) if needed. Without MPI, Gemini3D uses one CPU core only.
* [Ninja](https://ninja-build.org/) will build/rebuild much faster than GNU Make for any software project. `cmake -P scripts/install_ninja.cmake`

1. get the Gemini code

  ```sh
  git clone https://github.com/gemini3d/gemini3d.git

  cd gemini3d
  ```
2. Build Gemini and run self-test

  ```sh
  cmake -B build

  cmake --build build
  ```

Non-default [build options](./docs/Readme_cmake.md) may be used.

GEMINI has self tests that compare the output from a "known" test problem to a reference output.
To help ensure successful simulations, run the self-tests:

```sh
cd build

ctest
```

## How to setup a sim

1. make a [config.nml](./docs/Readme_input.md) with desired parameters for an equilibrium sim.
2. run the equilibrium sim:

    ```sh
    python -m gemini3d.run /path_to/config_eq.nml /path_to/sim_eq/
    ```
3. create a new config.nml for the actual simulation and run

    ```sh
    python -m gemini3d.run /path_to/config.nml /path_to/sim_out/
    ```

### Windows

Occasionally on Windows you may get a system error code `0xc0000005` when trying to run Gemini.
This typically requires rebooting the Windows computer.
If this is annoying, please let us know--it happens rarely enough that we're not sure if it's a Microsoft MPI bug or something else.

## Prerequisites

Gemini uses CMake build system to automatically build the entire software library stack,
checking for compatibility of pre-installed libraries such as Lapack, Scalapack and MUMPS.

### Libraries

Libraries are auto-built by Gemini when building gemini.bin.
These will generally yield faster Gemini runtime, since they were optimized for the CPU on your hardware.
If it's desired to use:

* system libraries: [PyGemini scripts/install_prereqs.py](https://github.com/gemini3d/pygemini)
* build/install libraries: `python -m gemini3d.prereqs`

## Known limitations and issues of GEMINI

1. Generating equilibrium conditions can be a bit tricky with curvilinear grids.  A low-res run can be done, but it will not necessary interpolate properly onto a finer grid due to some issue with the way the grids are made with ghost cells etc.  A workaround is to use a slightly narrower (x2) grid in the high-res run (quarter of a degree seems to work most of the time).
2. Magnetic field calculations on an open 2D grid do not appear completely consistent with model prototype results; although there are quite close.  This may have been related to sign errors in the FAC calculations - these tests should be retried at some point.
3. Occasionally MUMPS will throw an error because it underestimated the amount of memory needed for a solve.  If this happens a workaround is to add this line of code to the potential solver being used for your simulations.  If the problem persists try changing the number to 100.

    ```fortran
    mumps_par%ICNTL(14)=50
    ```
4. There are potentially some issues with the way the stability condition is evaluated, i.e. it is computed before the perp. drifts are solved so it is possible when using input data to overrun this especially if your target CFL number is &gt; 0.8 or so.  Some code has been added as of 8/20/2018 to throttle how much dt is allowed to change between time steps and this seems to completely fix this issue, but theoretically it could still happen; however this is probably very unlikely.
5. Occasionally one will see edge artifacts in either the field -aligned currents or other parameters for non-periodic in x3 solves.  This may be related to the divergence calculations needed for the parallel current (under EFL formulation) and for compression calculations in the multifluid module, but this needs to be investigated further...  This do not appear to affect solutions in the interior of the grid domain and can probably be safely ignored if your region of interest is sufficiently far from the boundary (which is always good practice anyway).

## Command-line options

By default, only the current simulation time and a few other messages are shown to keep logs uncluttered.
gemini.bin command line options include:

`-d` | `-debug`
: print verbosely -- could be 100s of megabytes of text on long simulation for advanced debugging.

`-nooutput`
: do not write data to disk. This is for benchmarking file output time, as the simulation output is lost, so this option would rarely be used.

`-out_format`
: normally Gemini reads and writes data files in the same format (HDF5, NetCDF4). This option allow one to read in one format (say NetCDF4) while writing HDF5.

* `h5`: HDF5 output (most commonly used)
* `nc`: NetCDF4 output
* `dat`: raw binary output (not recommended, doesn't support newer features)

`-manual_grid <# x2 images> <# x3 images>`
: forces the code to adopt a specific domain decomposition in x2 and x3 by using the integers given.  If not specified the code will attempt to find its own x2,x3 decomposition.  The number of grid points in x2 and x3 must be evenly divisible by the number of user-specified images in each direction, respectively.

`-dryrun`
: only run the first time step, do not write any files. This can be useful to diagnose issues not seen in unit tests, particularly issues with gridding. It runs in a few seconds or less than a minute for larger sims, something that can be done before queuing an HPC job.

### Number of MPI processes

In general for MPI programs and associated simulations, there may be a minimum number of MPI processes and/or integer multiples that must be met.
The build system generation process automatically sets the maximum number of processes possible based on your CPU core count and grid size.

This can also be done via `python -m gemini3d.run -np` options.

```sh
mpiexec -np <number of processors>  build/gemini.bin <output directory>
```

for example:

```sh
mpiexec -np 4 build/gemini.bin ~/mysim3d/arecibo
```

## Input file format

See [Readme_input](./docs/Readme_input.md)

## Loading and plotting output

GEMINI uses Python for essential interfaces, plotting and analysis.
Matlab scripts relevant to Gemini to
[mat_gemini repo](https://github.com/gemini3d/mat_gemini).

Only the essential scripts needed to setup a simple example, and plot the results are included in the main GEMINI repository.
The [Gemini-scripts](https://github.com/gemini3d/GEMINI-scripts)
and
[Gemini-examples](https://github.com/gemini3d/GEMINI-examples)
contain scripts used for various published and ongoing analyses.

See [Readme_output](./docs/Readme_output.md) for a description of how to load the simulation output files and the different variable names, meanings, and units.

An auxiliary program, magcalc.f90, can be used to compute magnetic field perturbations from a complete disturbance simulation.  See [Readme_magcalc](./docs/Readme_magcalc.md) for a full description of how this program works.


## List of other associated Readmes

1. [Readme_output](./docs/Readme_output.md)
2. [Readme_input](./docs/Readme_input.md)
3. [Readme_compilers](./docs/Readme_compilers.md)
4. [Readme_cmake](./docs/Readme_cmake.md)
5. [Readme_docs](./docs/Readme_docs.md)
6. [Readme_mpi](./docs/Readme_mpi.md)
7. [Readme_magcalc](./docs/Readme_magcalc.md)
8. [Readme_VEGA](./docs/Readme_VEGA.md)

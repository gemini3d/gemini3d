# GEMINI

[![DOI](https://zenodo.org/badge/146920930.svg)](https://zenodo.org/badge/latestdoi/146920930)
[![ci](https://github.com/gemini3d/gemini3d/actions/workflows/ci.yml/badge.svg)](https://github.com/gemini3d/gemini3d/actions/workflows/ci.yml)
[![oneapi-linux](https://github.com/gemini3d/gemini3d/actions/workflows/oneapi-linux.yml/badge.svg)](https://github.com/gemini3d/gemini3d/actions/workflows/oneapi-linux.yml)

The GEMINI model (*G*eospace *E*nvironment *M*odel of *I*on-*N*eutral *I*nteractions) is a three-dimensional ionospheric fluid-electrodynamic model written (mostly) in object-oriented fortran (2008+ standard).  GEMINI is used for various scientific studies including:

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

Generally, the Git `main` branch has the current development version and is the best place to start, while more thoroughly-tested releases happen regularly.
Specific releases corresponding to published results are generally noted in the corresponding journal article.

## Bug Reporting

The GEMINI development teams values input from our user community, particulary in the form of reporting of errors.  These allow us to insure that the code functions properly for a wider range of conditions, platforms, and use cases than we are otherwise able to directly test.

Please open a
[GitHub Issue](https://github.com/gemini3d/gemini/issues)
if you experience difficulty with GEMINI.  Try to provide as much detail as possible so we can try to reproduce your error.

## Platforms

Gemini is intended to be OS / CPU arch / platform / compiler agnostic.
Operating system support includes: Linux, MacOS, and Windows.
CPU arch support includes: Intel, AMD, ARM, IBM POWER, Cray and more.
GEMINI can run on hardware ranging from a Raspberry Pi to laptop to a high-performance computing (HPC) cluster.
Generally speaking one can run large 2D or modest resolution 3D simulations (less than 10 million grid points) on a quad-core workstation, with some patience.

For large 3D simulations (many tens-to-hundreds of millions of grid points), GEMINI is best run in a cluster environment or a very "large" multi-core workstation (e.g. 16 or more cores).
Runtime depends heavily on the grid spacing used, which determines the time step needed to insure stability,
For example we have found that a 20M grid point simulations takes about  4 hours on 72 Xeon E5 cores.  200M grid point simulations can take up to a week on 256 cores.
It has generally been found that acceptable performance requires > 1GB memory per core; moreover, a large amount of storage (hundreds of GB to several TB) is needed to store results from large simulations.

## Quick start

To build Gemini and run self-tests takes about 10 minutes on a laptop.
Gemini3D uses several external libraries that are built as a required one-time procedure.
Gemini3D works "offline" that is without internet once initially setup.

Requirements:

* C, C++ and Fortran compiler. See [compiler help](./docs/Readme_compilers.md) for optional further details.
  * GCC &ge; 9 with OpenMPI or MPICH
  * Clang with OpenMPI
  * Intel oneAPI
  * Cray with GCC or Intel oneAPI backend
* Python and/or MATLAB for scripting front- and back-ends
* CMake: if your CMake is too old, [download](https://cmake.org/download/) or `python -m pip install cmake`
* MPI: any of OpenMPI, IntelMPI, MPICH, MS-MPI. See [MPI help](./docs/Readme_mpi.md) if needed. Without MPI, Gemini3D uses one CPU core only, which runs much more slowly than with MPI.

### Gemini3D setup

Install Gemini3D prerequisite libraries. This is a one-time process used by any Gemini3D builds you do (or other programs). If your Python is too old, it will also install a local Python interpreter.

```sh
git clone https://github.com/gemini3d/external.git

cmake -P external/build-online.cmake
# installs under ~/libgem_gnu by default
```

---

Set environment variables `CMAKE_PREFIX_PATH` and edit PATH environment variable as follows.
On **Linux** add to ~/.bashrc, or on **MacOS** add to ~/.zshrc:

```sh
export CMAKE_PREFIX_PATH=~/libgem_gnu
```

---

Build the Gemini3D code

    ```sh
    git clone https://github.com/gemini3d/gemini3d.git

    cd gemini3d

    cmake -B build -DCMAKE_PREFIX_PATH=~/libgem

    cmake --build build --parallel
    ```

Non-default [build options](./docs/Readme_cmake.md) may be used.

GEMINI has self tests that compare the output from a "known" test problem to a reference output.
To verify your GEMINI build, run the self-tests.

```sh
ctest --test-dir build
```

### Offline HPC batch CTest

Note: some HPC systems only have internet when on a login node, but cannot run MPI simulations on the login node.
Batch sessions, including interactive, may be offline.
To run CTest in such an environment, download the data once from the login node:

```sh
ctest --test-dir build --preset download
```

then from an interactive batch session, run the tests:

```sh
ctest --test-dir build --preset offline
```

## GEMINI Numerical Library Dependencies

For various numerical solutions Gemini relies on:

* LAPACK
* scalapack
* MUMPS

For file input/output we also use:

* hdf5
* h5fortran
* zlib

## Running GEMINI from a Shell Environment

For basic operations the GEMINI main program simply needs to be run from the command line with arguments corresponding to to the number of processes to be used for the simulation, the location where the input files are and where the output files are to be written:

```sh
mpiexec -np <number of processors>  build/gemini.bin <output directory>
```

for example:

```sh
mpiexec -np 4 build/gemini.bin ~/mysim3d/arecibo
```

GEMINI can also be run via scripting frontends, e.g. `python -m gemini3d.run -np` options.

### Advanced Command Line Options

By default, only the current simulation time and a few other messages are shown to keep logs uncluttered.
gemini.bin command line options include:

`-d` | `-debug`
: print verbosely -- could be 100s of megabytes of text on long simulation for advanced debugging.

`-nooutput`
: do not write data to disk. This is for benchmarking file output time, as the simulation output is lost, so this option would rarely be used.

`-manual_grid <# x2 images> <# x3 images>`
: forces the code to adopt a specific domain decomposition in x2 and x3 by using the integers given.  If not specified the code will attempt to find its own x2,x3 decomposition.  The number of grid points in x2 and x3 must be evenly divisible by the number of user-specified images in each direction, respectively.

`-dryrun`
: only run the first time step, do not write any files. This can be useful to diagnose issues not seen in unit tests, particularly issues with gridding. It runs in a few seconds or less than a minute for larger sims, something that can be done before queuing an HPC job.


## Running GEMINI through Scripting Environments

If you prefer to issue the GEMINI run command through a scripting environment you may do so (via python) in the following way:

1. make a [config.nml](./docs/Readme_input.md) with desired parameters for an equilibrium sim.
2. run the equilibrium sim:

    ```sh
    python -m gemini3d.run /path_to/config_eq.nml /path_to/sim_eq/
    ```
3. create a new config.nml for the actual simulation and run

    ```sh
    python -m gemini3d.run /path_to/config.nml /path_to/sim_out/
    ```

## Input file format

See [Readme_input](./docs/Readme_input.md) for details on how to prepare input data for GEMINI.  Generally speaking there are python and MATLAB scripts available in the mat_gemini and pygemini repositories that will save data in the appropriate format once generated.

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


## Computing Magnetic Field Perturbations

An auxiliary program, magcalc.f90, can be used to compute magnetic field perturbations from a complete disturbance simulation.  See [Readme_magcalc](./docs/Readme_magcalc.md) for a full description of how this program works.


## List of other associated Readmes

1. [Readme_output](./docs/Readme_output.md) - information about data included in the output files of a GEMINI simulation
2. [Readme_input](./docs/Readme_input.md) - information on how input files should be prepared and formatted.
3. [Readme_compilers](./docs/Readme_compilers.md) - details regarding various compilers
4. [Readme_cmake](./docs/Readme_cmake.md) - cmake build options
5. [Readme_docs](./docs/Readme_docs.md) - information about model documentation
6. [Readme_mpi](./docs/Readme_mpi.md) - help with mpi-related issues
7. [Readme_magcalc](./docs/Readme_magcalc.md) - some documentation for the magnetic field calculation program
8. [Readme_VEGA](./docs/Readme_VEGA.md) - information on how to deploy and run GEMINI on ERAU's VEGA HPC system.
9. [Readme_prereqs](./docs/Readme_prereqs.md) - details on how to install prerequisites on common platforms.


## Known limitations and issues of GEMINI

1. Generating equilibrium conditions can be a bit tricky with curvilinear grids.  A low-res run can be done, but it will not necessary interpolate properly onto a finer grid due to some issue with the way the grids are made with ghost cells etc.  A workaround is to use a slightly narrower (x2) grid in the high-res run (quarter of a degree seems to work most of the time).
2. Magnetic field calculations on an open 2D grid do not appear completely consistent with model prototype results; although there are quite close.  This may have been related to sign errors in the FAC calculations - these tests should be retried at some point.
3. Occasionally MUMPS will throw an error because it underestimated the amount of memory needed for a solve.  If this happens a workaround is to add this line of code to the potential solver being used for your simulations.  If the problem persists try changing the number to 100.

    ```fortran
    mumps_par%ICNTL(14)=50
    ```
4. There are potentially some issues with the way the stability condition is evaluated, i.e. it is computed before the perp. drifts are solved so it is possible when using input data to overrun this especially if your target CFL number is &gt; 0.8 or so.  Some code has been added as of 8/20/2018 to throttle how much dt is allowed to change between time steps and this seems to completely fix this issue, but theoretically it could still happen; however this is probably very unlikely.
5. Occasionally one will see edge artifacts in either the field -aligned currents or other parameters for non-periodic in x3 solves.  This may be related to the divergence calculations needed for the parallel current (under EFL formulation) and for compression calculations in the multifluid module, but this needs to be investigated further...  This do not appear to affect solutions in the interior of the grid domain and can probably be safely ignored if your region of interest is sufficiently far from the boundary (which is always good practice anyway).
6. Occasionally on Windows you may get a system error code `0xc0000005` when trying to run Gemini.
This typically requires rebooting the Windows computer.
If this is annoying, please let us know--it happens rarely enough that we're not sure if it's a Microsoft MPI bug or something else.

# Installing GEMINI Prerequisites

As indicated in the main readme there are a number of prerequisites for GEMINI.  This readme covers some details of how to install these on various platforms.  We focus here on ease of use and describe the "least effort" methods for getting needed software.

## Verified local environment

The repository provides `scripts/install-local.sh` for Debian/Ubuntu, macOS and
Ubuntu inside Windows WSL. Run it from any directory using its absolute path.
Pass `--system-deps` to explicitly authorize native package installation with
APT or Homebrew. Without that flag it uses your existing system/HPC toolchain.
Other Linux distributions require equivalent native packages installed manually.

The installer creates an isolated Python environment, build tree and installation
under `build/local/` (override with `--root /absolute/path`). It installs the pinned
packages in `scripts/requirements-local.txt`, resolves the existing hash-pinned
native dependencies through CMake, builds the model, runs the standalone native
HDF5/MPI/MUMPS checks and unit tests, installs it, and checks the installed executable. Nothing is
installed into the system Python. A failed step returns a nonzero status and
leaves no successful `environment.json` record.

Python 3.11 or newer is required. On macOS, install Xcode command-line tools and
Homebrew first; `--system-deps` supplies Python 3.12 and the native libraries.
On Debian/Ubuntu, the distribution must supply Python 3.11+ (for example Debian 12
or Ubuntu 24.04); an older system Python must be upgraded explicitly.
`PYTHON=/absolute/path/to/python3` selects the bootstrap interpreter.
For Windows, `scripts/install-local.ps1 -SystemDeps` invokes the same installer
in the installed Ubuntu WSL distribution; `-Distribution` selects another WSL
distribution. It does not provision WSL or claim native Windows qualification.
The PowerShell entry point accepts `-Root`, `-Jobs`, `-BuildType Debug`,
`-ReferenceTests` and `-SourceCache` as equivalents of the shell options below.
`-Root` and `-SourceCache` are paths inside the selected WSL distribution, not
Windows drive paths. The Windows installation workflow invokes this PowerShell
entry point rather than bypassing it.

`--build-type Debug` enables the existing Debug checks. `--jobs` limits build/test
parallelism. Add `--reference-tests` to run **all registered tests**, including
reference downloads and simulations; the default unit-only installation is
**not** release or scientific qualification. Initial setup needs network access.
HPC/offline users should follow the existing offline dependency/reference-cache
instructions and provide their own compatible MPI/compiler environment.
`--source-cache /absolute/cache/sources.cmake` uses the hash-verified native
source cache produced by `scripts/offline_libraries.cmake`. Python wheels must
also be available to pip (for example through an administrator-provisioned local
index); this flag alone does not make an installation fully offline.

The successful `build/local/environment.json` records the source revision and
dirty status, per-file source hashes, native dependency pins/configuration,
Python packages, test inventory, verification scope and executable hash.
Inventory is not a claim that every file or scientific mode has been reviewed.
Enabled MSIS2/HWM14 parameter files are also required and hashed in the installed
resource directory, and their hashes are checked before running.

Use `python3 /absolute/path/to/scripts/local_environment.py check` to recheck an
installation, or `run --case /absolute/path/to/case --ranks 2` to launch it.
Pass native solver options after `--`, such as `-- -dryrun`. The runner checks the
recorded executable, Python package environment and MPI availability first; it
does not clean simulation output, download dependencies, or request privileges.
Missing or changed dependencies require an explicit installer rerun. Use the
same `--root` for installation, checks and runs. Direct native executable launches
remain supported but do not invoke this Python environment check.
The runner uses the MPI launcher selected by CMake and rejects a changed launcher
or dynamic-library/MPI-module environment rather than mixing MPI implementations.

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

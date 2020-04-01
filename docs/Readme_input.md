# Input configuration

Gemini uses Fortran 95 standard NAMELIST files.
Gemini expects that the files will be named like `inputs/config.nml`.

Each simulation needs an input file that specifies location of initial conditions and other pertinent information for the simulation.
Example config.nml are under the [tests/](./tests/) directory once you've built Gemini and run the self-tests.
Each subdirectory is a separate example usage of GEMINI for a particular problem.

A large number of examples (in addition to those included in the main repo) are included in the
[GEMINI-scripts](https://github.com/gemini3d/GEMINI-scripts)
repository.

## Running with different boundary and initial conditions:

GEMINI requires both initial and boundary conditions to run properly.  Specifically the user must provide a complete initial ionospheric state (density, drift, and temperature for all ionospheric species), along with boundary conditions for the electric potential (in 2D this are the top, bottom, and side potentials; in 3D the topside current density and side wave potentials).  Fluid state variables are given free-flow boundary conditions at the edges of the simulation grid.  The `io` module contains code dealing with input of initial state from file and the `potential_comm` and `potentialBCs_mumps` modules contains contains code dealing with boundary condition input.

There are presently two ways in which the boundary and initial conditions can be set for GEMINI:  subroutine-based input and file-based input.

PLEASE NOTE that future releases will use Fortran 2008 `submodule`, likely completely removing the option for subroutine-based initial and boundary conditions.

### Subroutine-based input (*not recommended* and to be deprecated in a future release):

There are two subroutines that can be modified by the user to provide boundary conditions to the code; these are described below. Note that, if any of these are changed, the code needs to be recompiled.

`./ionization/boundary_conditions/precipBCs_mod.f90` - the function `precipBCs' specifies the pattern of electron precipitation, including characteristic energy and total energy flux, over top of grid.  If the user does not specify an input file for precipitation boundary conditions in `config.ini`, then this subroutine will be called to set the boundary.

`./numerical/potential/boundary_conditions/potentialBCs_mumps.f90` - boundary conditions for the electric potential or field-aligned current.  The type of input that is being used is specified by the flags in the `config.ini` file for the simulation.  This subroutine will only be called if the user has not specified an input file containing boundary conditions.

By default these subroutines will be used for boundary conditions if file input is not specified in the config.ini input file.  The base GEMINI sets these to be zero potential (or current) and some negligible amount of precipitation.  Note that if you write over these subroutines then the code will use whatever you have put into them if file input is not specified.  This can lead to unintended behavior if ones modifies these and then forgets since the code will continue to use the modifications instead of some baseline.  Because of this issue, and the fact that GEMINI must be rebuilt every time these subroutines are changed, this method of boundary condition input is going to be removed.

### File-based input (*recommended*)

The file input is enabled by the appropriate flags (flagprecfile and flagE0file) set in the input `config.nml` file (see Section entitled "Input file format" above).
All examples included in `initialize/` in both the GEMINI and GEMINI-scripts repositories use this method for setting boundary conditions.
Note that the user can specify the boundary condition on a different grid from what the simulation is to be run with; in this case GEMINI will just interpolate the given boundary data onto the current simulation grid.

### Initial conditions

GEMINI needs density, drift, and temperature for each species that it simulations over the entire grid for which the simulation is being run as input.  Generally one will use the results of another GEMINI simulation that has been initialized in an arbitrary way but run for a full day to a proper ionospheric equilibrium as this input.  Any equilibrium simulation run this way must use full output (flagoutput=1 in the `config.ini`).  A useful approach for these equilibrium runs is to use a coarser grid so that the simulation completes quickly and then interpolate the results up to fine grid resolution.  An example of an equilibrium setup is given in `./initialize/2Dtest_eq`; note that this basically makes up an initial conditions (using `eqICs.m`) and runs until initial transients have settled.  An example of a script that interpolates the output of an equilibrium run to a finer grid is included with `./initialize/2Dtest`.

## Running one of the premade examples

Currently the main repo only includes the very basic 2Dtest and 3Dtest examples

## Creating a simulation

1. Create initial conditions for equilibrium simulation -  Several examples of equilibrium setups are included in the ./initialize directory; these end with `_eq`.  These are all based off of the general scripts `./setup/model_setup.m` and related scripts.  In general this equilbrium simulation will set the date, location, and geomagnetic conditions for the background ionospheric state for later perturbation simulations.
2. Run an equilibrium simulation at low resolution to obtain a background ionosphere.  See examples in ./initialize ending in `_eq`
3. Generate a grid - Several examples of grid generation scripts adapted to particular problems are given in the `initialize/` directory of the repo (see list above for an example).  In particular, for 2Dtest and 3Dtest there is a script that reads in an equilbirum simulation, creates a fine resolution mesh, and then interpolates the equilibrium data onto that fine mesh.
4. Interpolate the equilibrium results on to a high resolution grid and create new input files for full resolution - See examples in the ./initialize/ directories not ending in `_eq`.  These are all based off of the general `./setup/model_setup_interp.m` script.
5. Set up boundary conditions for potential, if required - see section of this document on boundary conditions
6. Set up precipitation boundary conditions, if required -  see section of this document on boundary conditions
7. Recompile the code with make *only if you are using subroutine based input and boundary conditions* (please note that this functionality will be removed in a later release).  If you are using file-based input then a rebuild is not necessary (this is another benefit of using file-based input)
8. Run your new simulation

## Running in two dimensions

The code determines 2D vs. 3D runs by the number of x2 or x3 grid points specified in the config.nml input file.
If the number of x2 grid points is 1, then a 2D run is executed (since message passing in the x3 direction will work normally).
If the number of x3 grid points is 1, the simulation will swap array dimensions and vector components between the x2 and x3 directions so that message passing parallelization still provides performance benefits.
The data will be swapped again before output so that the output files are structured normally and the user who is not modifying the source code need not concern themselves with this reordering.

## Example

```ini
&base
ymd = 2013,2,20               ! year, month, day
UTsec0 = 18000.0              ! UTsec0:  start time in UT seconds
tdur = 300.0                  ! tdur:  duration of simulation in seconds
dtout = 60.0                  ! dtout: how often to do file output
activ = 108.9, 111.0, 5       ! activ:  f107a,f107,Ap
tcfl = 0.9                    ! tcfl:  target cfl number
Teinf = 1500.0                ! Teinf:  exospheric electron temperature
/

&flags
potsolve = 1      ! solve electrodynamics:   0 - no; 1 - electrostatic; 2 - inductive
flagperiodic = 0
flagoutput = 1
flagcap = 0
flagdneu = 0
flagprecfile = 1
flagE0file = 1
flagglow = 1
/

&files
file_format = 'h5'
indat_size = 'tests/data/test3d_glow/inputs/simsize.h5'
indat_grid = 'tests/data/test3d_glow/inputs/simgrid.h5'
indat_file = 'tests/data/test3d_glow/inputs/initial_conditions.h5'
/

&setup
glat = 67.11
glon = 212.95
xdist = 200e3              ! eastward distance (meters)
ydist = 200e3               ! northward distance (meters)
alt_min = 80e3              ! minimum altitude (meters)
alt_max = 1000e3            ! maximum altitude (meters)
alt_scale = 13.75e3, 20e3, 200e3, 200e3  ! altitude grid scales (meters)
lxp = 20                    ! number of x-cells
lyp = 18                    ! number of y-cells
Bincl = 90                  ! geomagnetic inclination
nmf = 5e11
nme = 2e11
precip_latwidth = 0.25
precip_lonwidth = 0.25
Etarg = 50e-3   ! V/m
Efield_fracwidth = 0.142857
eqdir = 'tests/data/test3d_eq'
/

&precip
dtprec = 5.0
prec_dir = 'tests/data/test3d_glow/inputs/prec_inputs/'
/

&efield
dtE0 = 1.0
E0_dir = 'tests/data/test3d_glow/inputs/Efield_inputs/'
/

&glow
dtglow = 5.0
dtglowout = 60.0            ! dtglowout: how often to do Glow file output
/
```
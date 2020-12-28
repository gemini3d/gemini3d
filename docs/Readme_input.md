# GEMINI Input

In addition to command line options (see main README), GEMINI requires input file information to specify:

1. simulation configuration
2. grid size file
3. grid data strcuture file
4. neutral inputs
5. precipitation inputs
6. electric field inputs
7. initial conditions

## 1. Simulation configuration file

Gemini uses Fortran 95 standard NAMELIST files for the input configuration files.
Gemini will search the input directory location for files will be named like `inputs/config.nml` or `config.nml`.

Each simulation needs an input file that specifies location of initial conditions and other pertinent information for the simulation.
Example config.nml are under the [tests/](./tests/) directory once you've built Gemini and run the self-tests.
Each subdirectory is a separate example usage of GEMINI for a particular problem.

A large number of examples (in addition to those included in the main repo) are included in the
[GEMINI-examples](https://github.com/gemini3d/GEMINI-examples)
repository.

### Example config.nml input file

Note that most simulations will not use all of the input options shown here.  Source code reading in these parameters is in `config_nml.f90`.

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

! Flags controlling various aspects of solve and output behavior
&flags
potsolve = 1                  ! solve electrodynamics:   0 - no; 1 - electrostatic; 2 - inductive
flagperiodic = 0              ! whether to consider the x3 dimension periodic:  0 - no; 1 - yes
flagoutput = 1                ! what information to put in output files:  1 - all state variables; 2 - averaged plasma parameters; 3 - electron density only
/

! Inputs file locations and format
&files
file_format = 'h5'                                                 ! format of the input files
indat_size = 'tests/data/test3d_glow/inputs/simsize.h5'
indat_grid = 'tests/data/test3d_glow/inputs/simgrid.h5'
indat_file = 'tests/data/test3d_glow/inputs/initial_conditions.h5'
/

! This is only used by some matlab and python input scripts, the main fortran code ignores it
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

! (optional - default off) Include neutral atmospheric perturbation inputs from another model/dataset
&neutral_perturb
flagdneu = 1                       ! on/off for neutral perturbations:  0 - off; 1 - on
interptype = 3                     ! how to treat the input neutral data:  1 - axisymmetric; 2 - Cartesian; 3 - 3D Cartesian
sourcemlat = 44.9397d0             ! magnetic latitude of the source location
sourcemlon = 328.7981d0            ! magnetic longitude of the source location
dtneu = 6d0                        ! time step between neutral file inputs
drhon = 2d3                        ! neutral grid step size in the radial or y-direction (meridional)
dzn = 2d3                          ! neutral grid step, vertical
dxn = 2d3                          ! (only required if 3D) neutral grid step in x-direction (zonal)
source_dir = '../simulations/input/mooreOK_neutrals/'
/

! (optional - default off) Include disturbance precipitation based on file inputs
&precip
flagprecfile = 1                   ! use precipitaiton file input:  0 - no; 1 - yes
dtprec = 5.0                       ! time step between precipitation file inputs
prec_dir = 'tests/data/test3d_glow/inputs/prec_inputs/'
/

! (optional - default off) Include electric field boundary condition inputs from a file
&efield
flagE0file = 1                     ! use electric field bounary condition file input:  0 - no; 1 - yes
dtE0 = 1.0                         ! time step between electric field file inputs
E0_dir = 'tests/data/test3d_glow/inputs/Efield_inputs/'
/

! (optional - default off) Use glow to compute impact ionization, Cartesian grids only
&glow
flagglow = 1                ! use glow?  0 - no; 1 - yes
dtglow = 5.0                ! how often to recall GLOW to compute ionization
dtglowout = 60.0            ! ow often to do Glow file output
/

! (optional - default off) Controlling background precipitation characteristics
&precip_BG
PhiWBG=1e-5                   ! total energy flux (mW/m^2)
W0BG=3e3                      ! characteristic energy (eV)
/

! (optional - default off) Leading order electrodynamics
%capacitance
flagcap = 2                   ! whether to use ionospheric capacitance in the solves:  0 - no; 1 - ionospheric part; 2 - ionospheric+magnetospheric parts
magcap =  30.0                ! magnetospheric capacitance (Farads)
/

! (optional - defult off) whether or not to do milestone outputs; only works with h5 file output
&milestone
mcadence=10                   ! number of outputs per milestone
/

! (optional - default on) estimate parallel currents or not?
&Jpar
flagJpar=.false.
/

! (optional - default off) include gravitational drift terms in the drift equation and potential equation solutions
&gravdrift
flaggravdrift=.true.
/

! (optional - default off) control equatorial ionization anomaly
&EIA
flagEIA = .true.              ! toggles EIA calculation
v0equator = 10.0              ! equatorial peak drift value
/

! (optional - default value of 2)
&diffusion
diffsolvetype = 2             ! type of diffusion solver to use:  1 - backward Euler; 2 - TRBDF2
/

! (optional - off by default)
&lagrangian
flaglagrangian=.true.         ! whether or not to have the grid drift at the ExB speed
/
```

## 2,3. Grid input files

One of the most complicated parts of setting up a new simulation is creating a grid.
Grids are generated from scripts external to the main fortran code and then passes into GEMINI as files.
Generally when setting up a grid, it is likely easiest to work from an existing example in
[GEMINI-examples](https://github.com/gemini3d/GEMINI-examples).
In the event that none of the examples suffice as a starting point, the details of grid creations are documented below.

Document grid creation details here...

Grid structures, once created, are written to a file using the matlab `writegrid.m` API to insure that they have the correct file structure and arrangement.  I.e.

```matlab
gemini3d.write.grid(params,xg)
```
where the `xg` variable is a structure containing all the expected grid elements (see below).  The parameters argument is a structure containing this following fields:

1.  params.indat\_size
2.  params.indat\_grid


The writegrid API creates a file with the grid data structure in it as well as a small file containing the size information.

### Grid structure requirements

Grid structures, variable `xg` in the example above shall have the following fields:  

```MATLAB 
"x1", "x1i", "dx1b", "dx1h", "x2", "x2i", "dx2b", "dx2h", "x3", "x3i", "dx3b", "dx3h", "h1", "h2", "h3", "h1x1i", "h2x1i", "h3x1i", "h1x2i", "h2x2i", "h3x2i", "h1x3i", "h2x3i", "h3x3i", "gx1", "gx2", "gx3", "Bmag", "I", "nullpts", "e1", "e2", "e3", "er", "etheta", "ephi", "r", "theta", "phi", "x", "y", "z"
```

### Visualizing the grid

Plotgrid...  But explain how to use it...

## Neutral data input files

The examples of specifying and saving input neutral data input files are provided in [GEMINI-scripts](https://github.com/gemini3d/gemini-scripts/tree/master/magic/). The folder contains examples for preparation of 2D, 2D-axisymmetric and 3D neutral input files.

### Neutral input data requirements

Neutral input file data shall contain neutral fluid velocities, volumetric perturbations in temperature, and number densities for [O], [N_2] and [O_2].
2D Cartesian neutral inputs should contain meridional and vertical fluid velocities; 2D axisymmetric neutral inputs should contain radial and vertical fluid velocities; for 3D GEMINI simulations - meridional, zonal and vertical fluid velocities.
Neutral particle temperature perturbations represent averaged values over all species.

## 5,6. Running with different boundary and initial conditions:

GEMINI requires both initial and boundary conditions to run properly.  Specifically the user must provide a complete initial ionospheric state (density, drift, and temperature for all ionospheric species), along with boundary conditions for the electric potential (in 2D this are the top, bottom, and side potentials; in 3D the topside current density and side wave potentials).
Fluid state variables are given free-flow boundary conditions at the edges of the simulation grid.
The `io` module contains code dealing with input of initial state from file and the `potential_comm` and `potentialBCs_mumps` modules contains contains code dealing with boundary condition input.

There are presently two ways in which the boundary and initial conditions can be set for GEMINI:  subroutine-based input and file-based input.

PLEASE NOTE that future releases will use Fortran 2008 `submodule`, likely completely removing the option for subroutine-based initial and boundary conditions.

### Subroutine-based input (*not recommended and to be deprecated in a future release*):

There are two subroutines that can be modified by the user to provide boundary conditions to the code; these are described below. Note that, if any of these are changed, the code needs to be recompiled.

`./ionization/boundary_conditions/precipBCs_mod.f90` - the function "precipBCs" specifies the pattern of electron precipitation, including characteristic energy and total energy flux, over top of grid.
If the user does not specify an input file for precipitation boundary conditions in `config.nml`, then this subroutine will be called to set the boundary.

`./numerical/potential/boundary_conditions/potentialBCs_mumps.f90` - boundary conditions for the electric potential or field-aligned current.  The type of input that is being used is specified by the flags in the `config.nml` file for the simulation.
This subroutine will only be called if the user has not specified an input file containing boundary conditions.

By default these subroutines will be used for boundary conditions if file input is not specified in the config.nml input file.
The base GEMINI sets these to be zero potential (or current) and some negligible amount of precipitation.
Note that if you write over these subroutines then the code will use whatever you have put into them if file input is not specified.
This can lead to unintended behavior if ones modifies these and then forgets since the code will continue to use the modifications instead of some baseline.
Because of this issue, and the fact that GEMINI must be rebuilt every time these subroutines are changed, this method of boundary condition input is going to be removed.

### File-based input (*recommended*)

The file input is enabled by the appropriate flags (flagprecfile and flagE0file) set in the input `config.nml` file (see Section entitled [Example input config.nml file](#Example-config.nml-input-file)).
All examples included in `initialize/` in both the GEMINI and GEMINI-scripts repositories use this method for setting boundary conditions.
Note that the user can specify the boundary condition on a different grid from what the simulation is to be run with; in this case GEMINI will just interpolate the given boundary data onto the current simulation grid.

### Electric field input files requirements

Electric field input files shall contain the following information:

### Precipitation input files requirements

Precipitation input files shall contain the following variables.
These variables are one element per grid cell of the inputs/precip/simgrid.h5 file.
Each time step file has these variables.

```
E0p: Energy bin
Qp: Particle flux in this energy bin
```

the inputs/precip/simgrid.h5 file contains these vector variables:

```
mlat: magnetic latitude  (-90, 90)
mlon: magnetic longitude (0, 360)
```

the inputs/precip/simsize.h5 file contains these scalar variables:

```
llat: number of latitude cells
llon: number of longitude cells
```

The number of cells in the precipitation files is in general different than the number of simulation cells.
The Fortran code interpolates the precipitation data in space and time.

## 7. Initial conditions

GEMINI needs density, drift, and temperature for each species that it simulations over the entire grid for which the simulation is being run as input.  Generally one will use the results of another GEMINI simulation that has been initialized in an arbitrary way but run for a full day to a proper ionospheric equilibrium as this input.  Any equilibrium simulation run this way must use full output (flagoutput=1 in the `config.nml`).  A useful approach for these equilibrium runs is to use a coarser grid so that the simulation completes quickly and then interpolate the results up to fine grid resolution.  An example of an equilibrium setup is given in `./initialize/2Dtest_eq`; note that this basically makes up an initial conditions (using `eqICs.m`) and runs until initial transients have settled.  An example of a script that interpolates the output of an equilibrium run to a finer grid is included with `./initialize/2Dtest`.

### Initial condition input file requirements

Initial condition input files shall contain:


## Suggested workflow for creating input file to run a simulation

1. Create initial conditions for equilibrium simulation -  Several examples of equilibrium setups are included in the ./initialize directory; these end with `_eq`.  These are all based off of the general scripts `./setup/model_setup.m` and related scripts.  In general this equilbrium simulation will set the date, location, and geomagnetic conditions for the background ionospheric state for later perturbation simulations.
2. Run an equilibrium simulation at low resolution to obtain a background ionosphere.  See examples in ./initialize ending in `_eq`
3. Generate a grid - Several examples of grid generation scripts adapted to particular problems are given in the `initialize/` directory of the repo (see list above for an example).  In particular, for 2Dtest and 3Dtest there is a script that reads in an equilbirum simulation, creates a fine resolution mesh, and then interpolates the equilibrium data onto that fine mesh.
4. Interpolate the equilibrium results on to a high resolution grid and create new input files for full resolution - See examples in the ./initialize/ directories not ending in `_eq`.  These are all based off of the general `./setup/model_setup_interp.m` script.
5. Set up boundary conditions for potential, if required - see section of this document on boundary conditions
6. Set up precipitation boundary conditions, if required -  see section of this document on boundary conditions
7. Recompile the code with make *only if you are using subroutine based input and boundary conditions* (please note that this functionality will be removed in a later release).  If you are using file-based input then a rebuild is not necessary (this is another benefit of using file-based input)
8. Run your new simulation

### A note on running two-dimensional simulations

The code determines 2D vs. 3D runs by the number of x2 or x3 grid points specified in the config.nml input file.
If the number of x2 grid points is 1, then a 2D run is executed (since message passing in the x3 direction will work normally).
If the number of x3 grid points is 1, the simulation will swap array dimensions and vector components between the x2 and x3 directions so that message passing parallelization still provides performance benefits.
The data will be swapped again before output so that the output files are structured normally and the user who is not modifying the source code need not concern themselves with this reordering.

# GEMINI Input

In addition to command line options (see main README), GEMINI requires input file information to specify:

1. [Simulation Configuration File](#sim_config_file)
2. [Grid Input Files](#grid_input_files)
3. [Neutral Inputs](#neutral_input_files)
4. [Boundary Conditions](#boundary_conditions)
5. [Initial Conditions](#initial_conditions)

<a name="sim_config_file"></a>
## 1. Simulation Configuration File

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

```nml
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
flagperiodic = 0              ! whether to consider the x3 dimension periodic:  0 - no; nonzero - yes; 1 - yes and force periodicity in glat/glon across x3 (good for instability simulations needed a uniform background neutral atmosphere and SZA)
flagoutput = 1                ! what information to put in output files:  1 - all state variables; 2 - averaged plasma parameters; 3 - electron density only
/

! Inputs file locations and format
&files
file_format = 'h5'                                                 ! format of the input files
indat_size = 'test_data/test3d_glow/inputs/simsize.h5'
indat_grid = 'test_data/test3d_glow/inputs/simgrid.h5'
indat_file = 'test_data/test3d_glow/inputs/initial_conditions.h5'
/

! This is only used by some matlab and python input scripts, the main fortran code ignores it
&setup
glat = 67.11
glon = 212.95
xdist = 200e3              ! eastward distance (meters)
ydist = 200e3               ! northward distance (meters)
alt_min = 80e3              ! minimum altitude (meters)
alt_max = 1000e3            ! maximum altitude (meters)
alt_scale = 13.75e3, 20e3, 200e3, 200e3  ! parameters controlling nonuniform x1 grid, these are Cartesian and thus altitude
!> alt_scale = 10e3, 8e3, 500e3, 150e3   ! this is a high-resolution x1 grid; this one has been the default most most of Zettergren's publications that use Cartesian coordinates
x2parms = 400e3,18.8e3,50e3,100e3        ! parameters controlling nonuniform x2 grid
x3parms = 400e3,1.625e3,18.5e3,50e3.     ! parameters controlling nonniform x3 grid
lxp = 20                    ! number of x-cells
lyp = 18                    ! number of y-cells
Bincl = 90                  ! geomagnetic inclination
nmf = 5e11
nme = 2e11
precip_latwidth = 0.25
precip_lonwidth = 0.25
Etarg = 50e-3   ! V/m
Efield_fracwidth = 0.142857
eqdir = 'test_data/test3d_eq'
/

! (optional - default off) Include neutral atmospheric perturbation inputs from another model/dataset
&neutral_perturb
flagdneu = 1                       ! on/off for neutral perturbations:  0 - off; 1 - on
interptype = 3                     ! how to treat the input neutral data:  0 - 2D Cartesian; 1 - 2D axisymmetric; 3 - 3D Cartesian geomagnetic; 4 - 3D Cartesian geographic
sourcemlat = 44.9397d0             ! magnetic latitude of the source location
sourcemlon = 328.7981d0            ! magnetic longitude of the source location
dtneu = 6d0                        ! time step between neutral file inputs
drhon = 2d3                        ! neutral grid step size in the radial or y-direction (meridional)
dzn = 2d3                          ! neutral grid step, vertical
dxn = 2d3                          ! (only required if 3D) neutral grid step in x-direction (zonal)
source_dir = '../simulations/input/mooreOK_neutrals/'
/

! (optional - default off) allow neutral background atmosphere to change during simulation
&neutral_BG
flagneuBG = .true.       ! on or off
dtneuBG = 1800.0         ! how often to call MSIS and HWM (if applicable)
msis_version = 0         ! 0 or 20; which MSIS version to use, MSIS00 or MSIS 2.0
/

! (optional - default off) Include disturbance precipitation based on file inputs
&precip
flagprecfile = 1                   ! use precipitaiton file input:  0 - no; 1 - yes
dtprec = 5.0                       ! time step between precipitation file inputs
prec_dir = 'test_data/test3d_glow/inputs/prec_inputs/'
/

! (optional - default off) Include electric field boundary condition inputs from a file
&efield
flagE0file = 1                     ! use electric field boundary condition file input:  0 - no; 1 - yes
dtE0 = 1.0                         ! time step between electric field file inputs
E0_dir = 'test_data/test3d_glow/inputs/Efield_inputs/'
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

! (optional - default off) whether or not to do milestone outputs; only works with h5 file output
&milestone
mcadence=10                   ! number of outputs per milestone
/

! (optional - default on)
&Jpar
flagJpar=.false.              ! estimate parallel currents or not?
/

! (optional - default off)
&gravdrift
flaggravdrift=.true.          ! include gravitational drift terms in the drift equation and potential equation solutions
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

! (optional - off by default)
&diamagnetic
flagdiamagnetic=.true.        ! whether or not to compute pressure terms in perp. momentum balance
/
```

<a name="grid_input_files"></a>
## 2. Grid Input Files

One of the most complicated parts of setting up a new simulation is creating a grid.
Grids are generated from scripts external to the main fortran code and then passes into GEMINI as files.
Generally when setting up a grid, it is likely easiest to work from an existing example in
[GEMINI-examples](https://github.com/gemini3d/GEMINI-examples).
In the event that none of the examples suffice as a starting point, the details of grid creations are documented below.

If using the ```gemini3d.model.setup``` interface in MATLAB or Python, grid creation is controlled through the config.nml file, namely the parameters ```alt_scale,x2parms,x3parms``` under the ```setup``` namelist.  Each of these parameters has a four elements as follows:

```nml
alt_scale =
!> Formula for grid step size:
!>     dalt = d(1) + d(2) * tanh((alt(end) - d(3)) / d(4));
<dzref>       ! reference dz
<A>           ! amp of tanh for degrading resolution, min step size is dzref - A, max is dzref + A
<z0>          ! distance from "bottom" of grid where we start to degrade resolution
<ell>         ! transition length of degradation - a tanh(z/ell) profile is used for dz
```

For the horizontal grid distances (note that below is *not* in proper nml format; use code above is copy/pasting):

```nml
x2,{3}parms =
!> Formula for grid step size:
!>   dx = dx0 + dxincr * (1/2+1/2*tanh((x(end)-x2)/ell));
<degdist>,    ! distance from boundary at which we start to degrade resolution
<dx0>,        ! min step size for grid
<dxincr>,     ! max step size *increase* for grid, max grid step size (at the edges) will be dx0 + dxincr
<ell>         ! transition length of degradation - a tanh(x/ell) profile is used for dx
```

If using the lower-level ```gemini3d.grid.cart3d``` utility, the parameters above are passed into this function via optional input structure fields (see source code for details).  Note that the hyperbolic tangent profiles in x2,3 are mirrored about the x2,3=0 point so that the grid spacing is symmetric about the origin.

Grid structures, once created, are written to a file using the matlab `writegrid.m` API to insure that they have the correct file structure and arrangement.  I.e.

```matlab
gemini3d.write.grid(params,xg)
```
where the `xg` variable is a structure containing all the expected grid elements (see below).  The parameters argument is a structure containing this following fields:

1.  params.indat\_size
2.  params.indat\_grid


The writegrid API creates a file with the grid data structure in it as well as a small file containing the size information.

### Grid structure/file requirements

Grid structures, variable `xg` in the example above shall have the following fields and sizes (```lx1``` is the length of the x1 grid, ```lx2``` is the length of the x2 grid, ```lx3``` is the length of the x3 grid).

```MATLAB
"lx"            ! (3) number of grid points (excluding ghost cells) in x1,2,3 directions, viz. [lx1,lx2,lx3]
"x1",           ! (lx1+4) x1 position variable, including ghost cells
"x1i",          ! (lx1+1) x1 cell interface positions for non-ghost cells
"dx1b",         ! (lx1+3) backward differences along the x1-coordinate, excluding first ghost cell
"dx1h",         ! (lx1) x1 midpoint differences for all non-ghost cells
"x2",           ! (lx2+4) x2 position variable, including ghost cells
"x2i",          ! (lx2+1) x2 cell interface positions for non-ghost cells
"dx2b",         ! (lx2+3) backward differences along the x2-coordinate, excluding first ghost cell
"dx2h",         ! (lx2) x2 midpoint differences for all non-ghost cells
"x3",           ! (lx3+4) x3 position variable, including ghost cells
"x3i",          ! (lx3+1) x3 cell interface positions for non-ghost cells
"dx3b",         ! (lx3+3) backward differences along the x3-coordinate, excluding first ghost cell
"dx3h",         ! (lx1) x1 midpoint differences for all non-ghost cells
"h1",           ! (lx1+4,lx2+4,lx3+4) metric factor for x1, including in ghost cells
"h2",           ! (lx1+4,lx2+4,lx3+4) metric factor for x2, including in ghost cells
"h3",           ! (lx1+4,lx2+4,lx3+4) metric factor for x3, including in ghost cells
"h1x1i",        ! (lx1+1,lx2,lx3) x1 metric factor at the x1 cell interfaces and x2,3 cell centers
"h2x1i",        ! (lx1+1,lx2,lx3) x2 metric factor at the x1 cell interfaces and x2,3 cell centers
"h3x1i",        ! (lx1+1,lx2,lx3) x3 metric factor at the x1 cell interfaces and x2,3 cell centers
"h1x2i",        ! (lx1,lx2+1,lx3) x1 metric factor at the x2 cell interfaces and x1,3 cell centers
"h2x2i",        ! (lx1,lx2+1,lx3) x2 metric factor at the x2 cell interfaces and x1,3 cell centers
"h3x2i",        ! (lx1,lx2+1,lx3) x3 metric factor at the x2 cell interfaces and x1,3 cell centers
"h1x3i",        ! (lx1,lx2,lx3+1) x1 metric factor at the x3 cell interfaces and x1,2 cell centers
"h2x3i",        ! (lx1,lx2,lx3+1) x2 metric factor at the x3 cell interfaces and x1,2 cell centers
"h3x3i",        ! (lx1,lx2,lx3+1) x3 metric factor at the x3 cell interfaces and x1,2 cell centers
"gx1",          ! (lx1,lx2,lx3) x1-component of the gravitational field over the grid
"gx2",          ! (lx1,lx2,lx3) x2-component of the gravitational field over the grid
"gx3",          ! (lx1,lx2,lx3) x3-component of the gravitational field over the grid
"Bmag",         ! (lx1,lx2,lx3) the magnitude of the magnetic field at grid points
"I",            ! (lx2,lx3) the inclination angle of the magnetic field vs. "horizontal" (x2,3) locations on the grid
"nullpts",      ! (lx1,lx2,lx3) a floating point flag (0 or 1) indicating whether the core model code should treat each location as null (not part of the computational domain) or valid (to be used in calculations).
"e1",           ! (lx1,lx2,lx3,3) unit vector vs. x1,2,3 in the x1 direction, components in ECEF dipole cartesian coordinates
"e2",           ! (lx1,lx2,lx3,3) unit vector in the x2 direction in ECEF dipole cartesian coordinates
"e3",           ! (lx1,lx2,lx3,3) unit vector in the x3 direction in ECEF dipole cartesian coordinates
"er",           ! (lx1,lx2,lx3,3) unit vector in the radial direction in ECEF dipole cartesian coordinates
"etheta",       ! (lx1,lx2,lx3,3) unit vector vs. x1,2,3 in the theta (zenith angle) direction, components in ECEF dipole cartesian coordinates
"ephi",         ! (lx1,lx2,lx3,3) unit vector vs. x1,2,3 in the phi (zonal) direction, components in ECEF dipole cartesian coordinates
"r",            ! (lx1,lx2,lx3) radial coordinate (ECEF) as a function of x1,2,3
"theta",        ! (lx1,lx2,lx3) theta (zenith angle) coordinate as a function of x1,2,3
"phi",          ! (lx1,lx2,lx3) phi coordinate as a function of x1,2,3
"x",            ! (lx1,lx2,lx3) x coordinate (ECEF) as a function of x1,2,3
"y",            ! (lx1,lx2,lx3) y coordinate (ECEF) as a function of x1,2,3
"z"             ! (lx1,lx2,lx3) z coordinate (ECEF) as a function of x1,2,3
"glat"          ! (lx1,lx2,lx3) geographic latitude as a function of x1,2,3
"glon"          ! (lx1,lx2,lx3) geographic longitude as a function of x1,2,3
"alt"           ! (lx1,lx2,lx3) altitude as a function of x1,2,3
```

Not all of these variables will be used by the core model code; nevertheless they should be included in grid files so that postprocessing code can function properly.  Due to the complicated nature of the grid structure fields, it is *highly* recommended that you use one of the existing functions or user interfaces to create a grid structure and file.  Note also that the number of variables being tracked by the grid means that the grid files will occupy a large amount of storage space, but this prevents the code(s) from having to recompute metric factors, etc.

### Visualizing the grid

Plotgrid...  But explain how to use it...

<a name="neutral_input_files"></a>
## 3. Neutral data input files

The examples of specifying and saving input neutral data input files are provided in [GEMINI-scripts](https://github.com/gemini3d/gemini-scripts/tree/master/magic/). The folder contains examples for preparation of 2D, 2D-axisymmetric and 3D neutral input files.

### Neutral input data requirements

Neutral input file data shall contain neutral fluid velocities, volumetric perturbations in temperature, and number densities for [O], [N<sub>2</sub>] and [O<sub>2</sub>].
2D Cartesian neutral inputs should contain meridional and vertical fluid velocities as:

```
'/dnOall'       % O perturbations
'/dnN2all'      % N2 perturbations
'/dnO2all'      % O2 perturbations
'/dvnrhoall'    % dvnrhoall - fluid velocity in meridional direction or radial in Axisymmetric simulations
'/dvnzall'      % dvnzall - fluid velocity in vertical direction
'/dTnall'       % Temperature perturbations
```
These are stored in files named with the datetime corresponding to each neutral input frame, i.e. ```YYYYMMDD_SSSSS.SSSSSS.h5``` ("S" represents UT seconds with microsecond accuracy).  These arrays as permuted/indexed as:

```
dnOall(altitude,horizontal)
```

A ```simsize.h5``` input file is also required and has two variables:

```
'/lx1'          % number of *horizontal* (either radial or meridional) neutral input grid points
'/lx2'          % number of *vertical* neutral input grid points
```

2D axisymmetric neutral inputs should contain radial and vertical fluid velocities, stored in the same variables as for the 2D Cartesian input above.

For 3D GEMINI simulations the velocities have meridional, zonal, and vertical components:

```
'/dn0all'         % O perturbations
'/dnN2all'        % N2 perturbations
'/dnO2all'        % O2 perturbations
'/dvnxall'        % Zonal fluid velocity (positive east)
'/dvnrhoall'      % Meridional fluid velocity (positive north)
'/dvnzall'        % Vertical velocity (positive upward)
'/dTnall'         % Temperature perturbations
```

These arrays are permuted as:

```
dnOall(altitude,zonal,meridional)
```

The 3D neutral input simsize file is organized as:

```
'/lx1'      % number of zonal grid points
'/lx2'      % number of meridional points
'/lx3'      % number of vertical grid points
```

Neutral particle temperature perturbations in these input files are averaged over all species.

<a name="boundary_conditions"></a>
## 4. Boundary Conditions

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
All examples included in `init/` in both the GEMINI and GEMINI-scripts repositories use this method for setting boundary conditions.
Note that the user can specify the boundary condition on a different grid from what the simulation is to be run with; in this case GEMINI will just interpolate the given boundary data onto the current simulation grid.

### Electric field input files requirement

Electric field boundary/background condition input files include the grid over which the user is specifying the precipitation (```simsize.h5```, and ```simgrid.h5```) and the individual frames of input precipitation parameters (named using the datetime of the frame, e.g. ```20110302_01800.000000.h5``` etc.).

The ```inputs/Efield/simsize.h5``` file contains these scalar variables:

```
"llat"        ! number of latitude cells
"llon"        ! number of longitude cells
```

The inputs/precip/simgrid.h5 file contains these vector variables:

```
"mlon"        ! (llon) magnetic longitude, units of degrees (0, 360)
"mlat".       ! (llat) magnetic latitude, units of degrees  (-90, 90)
```

Each frame electric field input file has variables:

```
"Exit"          ! (llon,llat) background electric field x-component (units of V/m in the x2 direction) vs. mlon,mlat
"Eyit"          ! (llon,llat) background electric field y-component (V/m in the x3 direction) vs. mlon,mlat
"Vmaxx1it"      ! (llon,llat) potential or FAC boundary condition at the location of maximum x1 vs. mlon,mlat
"Vmaxx2ist"     ! (llat) potential boundary condition at the location of maximum x2 vs. mlat
"Vmaxx3ist"     ! (llon) potential boundary condition at the location of maximum x3 vs. mlon
"Vminx1it"      ! (llon,llat) potential or FAC boundary condition at the location of minimum x1 vs. mlon,mlat
"Vminx2ist"     ! (llat) potential boundary condition at the location of minimum x2 vs. mlat
"Vminx3ist"     ! (llon) potential boundary condition at the location of minimum x3 vs. mlon
"flagdirich"    ! (1) whether to treat the data in Vmax{min}x1 arrays as potential (0 value) or FAC (1 value)
```

In cases where a quasi-electrodynamic solution is done both the electric field and potential will be assumed to go to zero at the lateral (x2,3) boundaries of the model - unless periodic x3 boundary conditions are chosen (these are set in the config.nml file).

### Precipitation input files requirement

Precipitation input files include the grid over which the user is specifying the precipitation (```simsize.h5``` and ```simgrid.h5```) and the individual frames of input precipitation parameters (named using the datetime of the frame, e.g. ```20110302_01800.000000.h5``` etc.).

The inputs/precip/simsize.h5 file contains these scalar variables:

```
"llat"        ! number of latitude cells
"llon"        ! number of longitude cells
```

The inputs/precip/simgrid.h5 file contains these vector variables:

```
"mlon"        ! (llon) magnetic longitude, units of degrees (0, 360)
"mlat"        ! (llat) magnetic latitude, units of degrees  (-90, 90)
```

These variables are one element per grid cell of the ```inputs/precip/simgrid.h5``` file.

Each frame precipitation input file has variables:

```
"E0p"         ! (llon,llat) Characteristic energy (Maxwellian differential number flux), units of eV
"Qp"          ! (llon,llat) Total energy flux in this energy bin, units of W/m^2
```

The number of cells in the precipitation files is in general different than the number of simulation cells.
The Fortran code interpolates the precipitation data in space and time.

<a name="initial_conditions"></a>
## 5. Initial conditions

GEMINI needs density, drift, and temperature for each species that it simulations over the entire grid for which the simulation is being run as input.  Generally one will use the results of another GEMINI simulation that has been initialized in an arbitrary way but run for a full day to a proper ionospheric equilibrium as this input.  Any equilibrium simulation run this way must use full output (flagoutput=1 in the `config.nml`).  A useful approach for these equilibrium runs is to use a coarser grid so that the simulation completes quickly and then interpolate the results up to fine grid resolution.  An example of a 2D equilibrium
[config.nml](https://github.com/gemini3d/gemci/tree/main/cfg/equilibrium/mini2dew_eq):
note that this basically makes up initial conditions and runs until initial transients have settled.
An example
[config.nml](https://github.com/gemini3d/gemci/tree/main/cfg/hourly/mini2dew_fang) interpolates the output of an equilibrium run to a finer grid.

### Initial condition input file requirements

Initial condition input files shall contain all input data needed to start a simulation including state variables for all plasma species (density, drift (parallel dimension), and temperature) in SI units.
These variables are to be organized as follows (```lsp``` is the number of species used in the simulations):

 ```
 "nsall"           ! (lx1,lx2,lx3,lsp) number density of each species over the grid
 "vs1all"          ! (lx1,lx2,lx3,lsp) drift velocity parallel to B over the grid (x1-direction)
 "Tsall"           ! (lx1,lx2,lx3,lsp) tmeperature for each species over the grid
 "Phiall"          ! (lx2,lx3) electric potential vs. x2 and x3 - may be omitted to default to zero
 ```

## Suggested workflow for creating input file to run a simulation

1. Create initial conditions for equilibrium simulation -  Several examples of equilibrium [config.nml](https://github.com/gemini3d/gemci/tree/main/cfg/equilibrium).  These are all based off of the general scripts `./setup/model_setup.m` and related scripts.  In general this equilbrium simulation will set the date, location, and geomagnetic conditions for the background ionospheric state for later perturbation simulations.
2. Run an equilibrium simulation at low resolution to obtain a background ionosphere.
3. Generate a grid - Several examples of grid generation scripts: see list above for an example.  From Matlab or Python, the gemini3d.model.setup() reads in an equilbirum simulation, creates a fine resolution mesh, and then interpolates the equilibrium data onto that fine mesh.
4. Interpolate the equilibrium results on to a high resolution grid and create new input files for full resolution.
5. Set up boundary conditions for potential, if required - see section of this document on boundary conditions
6. Set up precipitation boundary conditions, if required -  see section of this document on boundary conditions
7. Recompile the code with make *only if you are using subroutine based input and boundary conditions* (please note that this functionality will be removed in a later release).  If you are using file-based input then a rebuild is not necessary (this is another benefit of using file-based input)
8. Run your new simulation

## A note on running two-dimensional simulations

The code determines 2D vs. 3D runs by the number of x2 or x3 grid points specified in the config.nml input file.
If the number of x2 grid points is 1, then a 2D run is executed (since message passing in the x3 direction will work normally).
If the number of x3 grid points is 1, the simulation will swap array dimensions and vector components between the x2 and x3 directions so that message passing parallelization still provides performance benefits.
The data will be swapped again before output so that the output files are structured normally and the user who is not modifying the source code need not concern themselves with this reordering.

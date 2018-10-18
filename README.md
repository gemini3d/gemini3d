# GEMINI3D

The GEMINI model (*G*eospace *E*nvironment *M*odel of *I*on-*N*eutral *I*nteractions) is a three-dimensional ionospheric fluid-electrodynamic model used for various scientific studies including effects of auroras on the terrestrial ionosphere, natural hazard effects on the space environment, and effects of ionospheric fluid instabilities on radio propagation (see references section of this document for details).  The detailed mathematical formulation of GEMINI3D is included in `docs/`.  A subroutine-level set of documentation describing functions of individual program units is given via source code comments.  

The current version of the code uses generalized orthogonal curvilinear coordinates and has been tested with dipole and Cartesian coordinates.

We have prioritized ease of setup/install across a wide variety of computing systems.
Please open a [GitHub Issue](https://github.com/mattzett/gemini/issues) if you experience difficulty.

## Prerequisites

The CMake build system is the powerful and easy to use for large Fortran projects.
In general a recent CMake version is beneficial for Fortran builds.
CMake &ge; 3.11 is required, and easily installed without sudo on:

* Linux: use [cmake_setup.sh](https://github.com/scivision/cmake-utils)
* MacOS: `brew install cmake`
* [Windows](https://cmake.org/download)

### Compilers

MPI is woven throughout GEMINI3D so compiler wrappers `mpifort` or `mpiifort` can be often used.

* gfortran 4.8 - 7.3
* Intel `ifort`

### Libraries

used / tested versions include:

* OpenMPI 2.1 - 3.1
* MUMPS 5.1
* SCALAPACK 2.0
* LAPACK95 3.0  (optional)

### data file analysis

GEMINI3D `*.m` scripts require EITHER:

* GNU Octave &ge; 4.0
* Matlab &ge; R2007b

Generally, the Git `master` branch has the current development version and is the best place to start, while more thoroughly-tested releases happen occasionally.

## License

GEMINI3D is distributed under the Affero GNU public license (aGPL) version 3+.


## Suggested hardware

GEMINI3D should be run in a cluster environment or using a "large" multicore workstation, except in cases where 2D simulations are used. In general one could run large 2D or very small 3D simulations (not exceeding a few million grid points) on a quad-core workstation; anything larger than this needs to be run with 8-64 cores.  If resources and simulation geometry allow (note the parallelization is only along the 3rd dimensions (except for purely 2D runs), the optimal multi-node cluster setup will have about 6-7 GB memory per core.  A large amount of storage space is needed to store the results as large 3D simulations can generate 1-2 TB output.  Smaller 3D and 2D simulations can usually fit into tens of GB of storage space (note that the code does no compression on the output - to reduce the already significant output times).  


## Quick start
This method is tested on CentOS and Ubuntu.
This test runs a short demo, taking about 2-5 minutes on a typical Mac / Linux laptop, from scratch. 


1. get GEMINI code and install prereqs
   ```sh
   cd ~
   git clone https://github.com/mattzett/gemini
   cd gemini
   ```
2. Get the 2D test data: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1464915.svg)](https://doi.org/10.5281/zenodo.1464915)
3. compile and run GEMINI demo: 
   ```sh
   cd ~/gemini/objects

   cmake ..
   
   cmake --build .

   ctest -V
   ```
   
#### Build tips

* If CMake version too old, use [cmake_setup.sh](https://github.com/scivision/cmake-utils). This does NOT use `sudo`.
* If missing prereqs, try the `./install_prereqs.sh` script.
* If need to build libraries from source (e.g. because you don't have `sudo`) try `build_gfortran.sh` or `build_intel.sh` from the `fortran-libs` repo:
  ```sh
  cd ~
  git clone https://github.com/scivision/fortran-libs
  cd fortran-libs/LAPACK95
  make double -C SRC
  ```


### self-tests
GEMINI has self tests that compare the output from a "known" test problem to a reference output.  So running:
```sh
make test
```

1. executes 
   ```sh
   ./gemini initialize/2Dtest/config.ini /tmp/2d
   ```
2. uses GNU Octave (free lightweight Matlab clone) to compare with reference output using `tests/compare_all.m`:
   ```matlab
   compare_all(YourOutputDirectory, '../simulations/2Dtest_files/2Dtest_output')
   ```

### Ubuntu
Tested on Ubuntu 18.04 / 16.04.

If you have sudo (admin) access:
```sh
./install_prereqs.sh
```
Otherwise, ask your IT admin to install the libraries or 
[compile them yourself](https://github.com/scivision/fortran-libs) 
or consider Linuxbrew.


### CentOS
This is for CentOS 7, using "modules" for more recent libraries.
For the unavailable modules, 
[compile them yourself](https://github.com/scivision/fortran-libs)
```sh
module load git cmake mumps scalapack openmpi lapack metis

module load gcc/6.2.0
export CC=gcc CXX=g++ FC=gfortran
```

Try to compile gemini as above, then 
[build the missing libraries](https://github.com/scivision/fortran-libs).

Example:
```sh
cmake -DSCALAPACK_ROOT=~/fortran-libs/scalapack -DMETIS_ROOT=/share/pkg/metis/5.1.0/install ..
```

Alternatively for Intel Fortran:
```sh
FC=ifort cmake -DMETIS_ROOT=~/fortran-libs/metis ..
```

## Known limitations and issues of GEMINI

1)  Generating equilibrium conditions can be a bit tricky with curvilinear grids.  A low-res run can be done, but it will not necessary interpolate properly onto a finer grid due to some issue with the way the grids are made with ghost cells etc.  A workaround is to use a slightly narrower (x2) grid in the high-res run (quarter of a degree seems to work most of the time).
2)  Magnetic field calculations on an open 2D grid do not appear completely consistent with MATLAB model prototype results; although there are quite close.  This may have been related to sign errors in the FAC calculations - these tests should be retried at some point.  
3)  Occasionally MUMPS will throw an error because it underestimated the amount of memory needed for a solve.  If this happens a workaround is to uncomment (or add) this line of code to the potential solver being used for your simulations:
    ```fortran
    mumps_par%ICNTL(14)=50
    ```
If the problem persists try changing the number to 100. 
4)  There are potentially some issues with the way the stability condition is evaluated, i.e. it is computed before the perp. drifts are solved so it is possible when using input data to overrun this especially if your target CFL number is > 0.8 or so.  Some code has been added as of 8/20/2018 to throttle how much dt is allowed to change between time steps and this seems to completely fix this issue, but theoretically it may still happen.  
5)  Occassionally one will see edge artifacts in either the field -aligned currents or other parameters for non-periodic in x3 solves.  This may be related to the divergence calculations needed for the parallel current (under EFL formulation) and for compression calculations in the multifluid module, but this needs to be investigated further...


## Development roadmap

### To do list

### Future Code refactoring

Code duplication:
* Fair bit of code repetition in top-level electric field and precipitation interpolation routines
* Code duplication in electrodynamics module (haloing part should be written as a subroutine)
* Cleanup of BCs interpolation source files
* Axisymmetric and Cartesian interpolations should be combined (much code-sharing)

Research needed:  
* What's up with unallocated filenamefull in output\_root\_mpi???  f2003 feature that autoallocates strings???

Other development-related comments:
* There may be a performance boost by using the Fortran 2008 `contiguous` attribute on the `pointer` arrays where right now it is manually repacked--`contiguous` means we DON'T repack manually, the compiler will repack IF and ONLY IF it needs too.  We may get a performance boost by eliminating manual repacking and using `real, contiguous, pointer` instead. [Reference](https://modelingguru.nasa.gov/servlet/JiveServlet/previewBody/1527-102-1-2631/N1729-4.pdf) page 7.

### Feature requests
* (SOMEWHAT IMPORTANT) Parallel domain decomposition in x2 *and* x3 - this is a big task that is likely to be left aside until I can renew funding.  It's also questionable how useful it is at this point where my typical runs are 32-256 cores (although undoubtedly it may become useful for runs with thousands of cores).  I've found good speedup even dividing the x3 dimension into slabs 2 grid points wide; although that means passing essentially all the grid data around via mpi, the large number of operations per slab means that the effective overhead here is not too much to prevent this from being useful.  
* GLOW functionality merged into curvilinear branch
* (SOMEWHAT IMPORTANT) Periodically updating background neutral atmosphere - should really be done for simulations more than a few hours long but will affect performance
* (EFFICIENCY) Exclusion of null points from field aligned advection, thermal conduction, and source terms - could improve performance
* Possibly merge in P. Inchin's EIA changes (with appropriate flags)
* HDF5 file input and output
* Option to run the code in a single precision mode - would help with memory limited systems although it's not clear how this would impact numerics (I've never tested my methods in single precision)

### Plans for adding physics:

These are projects in progress involved GEMINI, you are encouraged to email M. Zettergren for more info if you have interest in using or collaborating on these so that we can efficiently combine efforts and avoid duplicative work.

* Resolved potential solutions - decimate parallel grid down to Farley mapping scale for perp resolution then so the solve on that coarse grid then interpolate back up to original grid.  I've had luck with MUMPS solves in reasonable time up to 300 x 300 x 15 grid points which is probably enough to do something interesting with appropriate periodic and lagrangian grids (moving at E x B).  
* Diamagnetic drift and perpendicular ambipolar fields - necessary for the smallest scales, e.g. less than 100 m
* Inclusion of suprathermal electron transport model for better specification of currents and ionization rates (G. Grubbs)


## Standard and style

GEMINI3D is Fortran 2018 compliant and uses two-space indents throughout (to accommodate the many, deeply nested loops).  If you modify the source code and plan to push back to the repository please *do not* use tabs.



## To build and run GEMINI3D:

    make clean
    make all
    mpirun -np <number of processors>  ./gemini_mumps <input config file> <output directory>

for example:  

    mpirun -np 4 ./gemini_mumps ./initialize/2Dtest/config.ini ~/simulations/2Dtest/

Note that the output *base* directory must already exist (‘simulations’ in previous example).  The source code consists of about ten module source files encapsulating various functionalities used in the model.  A diagram all of the modules and their function is shown in figure 1; a list of module dependencies can also be found in the Makefile source.

![Figure 1](doc/figure1.png)

<!-- ![Figure 2](doc/figure2.png) -->


## Verifying GEMINI build

Assuming you have built by
```sh
cd objects
cmake ..
cmake --build . 
```
you can run the self test with
```sh
ctest -V
```

## Input file format

Each simulation needs an input file that specifies location of initial conditions and other pertinent information for the simulation.  Numerous examples of these are included in the ./initialize directory; each subdirectory is a separate example usage of GEMINI for a particular problem.  The basic template for an input file (config.ini) file follows (please note that most use cases will not have all options activated as this example does).  
```
16,9,2015                             !dmy:  day,month,year
82473.0                               !UTsec0:  start time, UT seconds
1800.0                                !tdur:  duration of simulation in seconds
15.0                                  !dtout: how often (s) to do output
109.0,109.0,5.0                       !activ:  f107a,f107,Ap (81 day averaged f10.7, daily f10.7, and average Ap index)
0.9                                   !tcfl:  target cfl number (dimensionless - must be < 1d0 to insure stability)
1500.0                                !Teinf:  exospheric electron temperature, K (only used in open-grid simulations)
0                             	  !potsolve:  are we solving potential? (0=no; 1=-yes)
0                                     !flagperiodic:  do we interpret the grid as being periodic in the x3-direction?  (0=no; 1=yes)
2                                     !flagoutput:  what type of output do we do?  (2=ISR-like species-averaged plasma parameters; 3=electron density only; anything else nonzero=full output)
0                                     !flagcapacitance:  include inertial capacitance? (0=no; 1=yes; 2=yes+m'spheric contribution)
../zettergmdata/simulations/input/chile20153D_0.5_medhighres/chile20153D_0.5_medhighres_simsize.dat
../zettergmdata/simulations/input/chile20153D_0.5_medhighres/chile20153D_0.5_medhighres_simgrid.dat
../zettergmdata/simulations/input/chile20153D_0.5_medhighres/chile20153D_0.5_medhighres_ICs.dat
1                                     !are we applying neutral perturbations? (0=no; 1=yes).  If 0, the next five entries are skipped while reading this input file
1                                     !how doe we interpret the input neutral file geometry?  (0=Cartesian; anything else=axisymmetric)
-20.5706d0,359.4048d0                 !source mlat,mlon of disturbance (degrees magnetic lat,lon)
4d0                                   !time step between neutral input files
2d3,2d3                               !spatial resolutions in radial and vertical directions
../zettergmdata/simulations/chile2015_neutrals/
1                                     !flagprecfileinput:  for precipitation file input (0=no; 1=yes).  If 0, then next two entries of input file are skipped
1.0                                   !dtprec:  time (s) between precipitation input files
../simulations/isinglass_precipitation/
1                                     !flagE0fileinput:  flag for electric field file input (0-no; 1=yes).  If 0, next two entries of input file are skipped
10.0                                  !dtE0:  time (s) between electric field input files
../simulations/isinglass_fields/
```


## Running with different boundary and initial conditions:  

GEMINI requires both initial and boundary conditions to run properly.  Specifically the user must provide a complete initial ionospheric state (density, drift, and temperature for all ionospheric species), along with boundary conditions for the electric potential (in 2D this are the top, bottom, and side potentials; in 3D the topside current density and side wave potentials).  Fluid state variables are given free-flow boundary conditions at the edges of the simulation grid.  The `io` module contains code dealing with input of initial state from file and the `potential_comm` and `potentialBCs_mumps` modules contains contains code dealing with boundary condition input.  

There are presently two ways in which the boundary and initial conditions can be set for GEMINI:  subroutine-based input and file-based input.  Future releases will likely completely remove the option for subroutine-based initial and boundary conditions.  

### Subroutine-based input (*not recommended* and soon to be deprecated):  

There are two subroutines that can be modified by the user to provide boundary conditions to the code; these are described below. Note that, if any of these are changed, the code needs to be recompiled.  

./ionization/boundary\_conditions/precipBCs\_mod.f90 - the function `precipBCs' specifies the pattern of electron precipitation, including characteristic energy and total energy flux, over top of grid. 

./numerical/potential/boundary_conditions/potentialBCs_mumps.f90 - boundary conditions for the electric potential or field-aligned current.  The type of input that is being used is specified by the flags in the config.ini file for the simulation.  

By default these subroutines will be used for boundary conditions if file input is not specified in the config.ini input file.  These are initially set in the source code to be zero potential (or current) and some negligible amount of precipitation.  Note that if you write over these subroutines then the code will use whatever you have put into them if file input is not specified.  This can lead to unintended behavior if ones modifies these and then forgets since the code will continue to use the modifications instead of some baseline.  Because of this issue, the mode of specifying boundary conditions should probably be entirely removed in a later version of the code but for now it is being kept to maintain compatibility with some older projects.  

### File-based input (*recommended*)

An alternative is to use the file input option, which needs to be set up using MATLAB (or other) scripts.  To enable this type of input, the appropriate flags (flagprecfileinput and flagE0fileinput) need to be set in the input `config.ini` file (see Section entitled "Input file format" above).  Several examples of using file-based input are included in `initialize/`; e.g. see the `tohoku2011`, `GDI`, `KHI` and `isinglass` examples.  


## Running one of the premade examples:

Several different examples are included with the source code; although initial conditions for each must be generated by the user by running a corresponding equilibrium simulations which generates balanced initial conditions for a given date, time, etc.  These equilibrium runs generally are started with a made-up initial condition so there is a lot of initial settling before a sensible ionospheric state is achieved.  To account for this one usually needs to run for about 24 hours of simulation time to insure a set of state parameters that are a good representation of the ionosphere..  Each of these examples has its own initial and boundary conditions generation scripts which are stored in the appropriately named directories in the ./initialize/ directory, along with a config.ini file need as input to the simulation.  The generation scripts must be run in order to produce input grids and initial conditions for each simulation.

The examples are labeled:

* ARCS_eq - an equilibrium (eq) simulation that generates initial conditions (ICs) for the ARCS simulation described below.
* ISINGLASS_eq - an equilibrium simulation generating ICs for the ISINGLASS simulation
* tohoku20113D_eq - eq simulation for 3D Tohoku earthquake simulations
* tohoku20112D_eq - eq simulation for 2D Tohoku earthquake simulations
* RISR_eq - eq simulation for the GDI and KHI examples described below (location:  Resolute Bay ISR)
* nepal20152D_eq - eq simulation for the 2D nepal earthquake simulation
* 2Dtest_eq - eq simulation for 2D test case that can be run quickly to verify the code is working correctly.
* tohoku20113D_medres - A medium resolution simulation in 3D of the 2011 Tohoku earthquake.
* nepal20152D_highres - A high resolution simulation in 2D for the 2015 Nepal eartquake.
* tohoku20112D_highres - A 2D, high resolution simulation for the 2011 Tohoku earthquake.
* ARCS - 
* ISINGLASS - an example using model input derived from observations from the ISINGLASS sounding rocket campaign.  This example illustrates how to drive the model with data inputs.
* GDI_periodic_medres_fileinput - a simulation of gradient-drift instability illustrating the use of a periodic mesh
* KHI_periodic_highres_fileinput - a simulation of Kelvin-Helmholtz instability illustration periodic meshes and use of polarization current solver
* 2Dtest - a quick to run (3 minutes on quadcore) example for testing whether the code compiled and runs properly.

A fair bit of testing has been done on these, but there could still be problems so contact a developer if you are having issues with the examples.


## Creating a simulation

1)  Generate a grid - Several examples of grid generation scripts adapted to particular problems are given in the ./initialize directory of the repo (see list above for an example).  These are all based off of the general scripts:  ./setup/gridsplot.m OR ./setup/gridsplot\_map.m
2)  Create initial conditions for equilibrium simulation -  Several examples of equilibrium setups are included in the ./initialize directory; these end with "\_eq".  These are all based off of the general scripts ./setup/model\_setup.m and related scripts.  
3)  Run an equilibrium simulation at low resolution to obtain a background ionosphere.  See examples in ./initialize ending in "_eq"
4)  Interpolate the equilibrium results on to a high resolution grid and create new input files for full resolution - See examples in the ./initialize/ directories not ending in "\_eq"  These are all based off of the general ./setup/model\_setup\_interp.m script. 
5)  Set up boundary conditions for potential, if required - see section of this document on boundary conditions
6)  Set up precipitation boundary conditions, if required -  see section of this document on boundary conditions
7)  Recompile the code with make *only if you are using subroutine based input and boundary conditions* (please note that this functionality will be removed in a later release).  If you are using file-based input then a rebuild is not necessary (this is another benefit of using file-based input)
8)  Run your new simulation


## Running in two dimensions

The code determines 2D vs. 3D runs by the number of x2 or x3 grid points specified in the 'config.ini' input file.  If the number of x2 grid points is 1, then a 2D run is executed (since message passing in the x3 direction will work normally).  If the number of x3 grid points is 1, the simulation will swap array dimensions and vector components between the x2 and x3 directions so that message passing parallelization still provides performance benefits.  The data will be swapped again before output so that the output files are structured normally and the user who is not modifying the source code need not concern themselves with this reordering.


## Loading and plotting output

Either MATLAB or GNU/octave is required to load the output file via scripts in the ./vis directory (these scripts generally work on both 2D and 3D simulation results).  The results for an entire simulation can be plotted with 'plotall.m' (see source code for details), which also illustrates how to read in a sequence of files from a simulation.  This script prints a copy of the output plots into the simulation output directory.  Finer-level output control can be achieve by using the 'plotframe.m' and 'loadframe.m' scripts to plot and load data from individual simulation output frames, respectively.  

Output frames (corresponding to a snapshot of the ionospheric state at a particular time) are each stored in a different binary file named according to the date and time to which the frame corresponds; these may be read into the MATLAB GNU/octave workspace using the `loadframe.m` and `plotall.m` MATLAB functions.  There is also a `plotframe.m` function provided which will produce plots for a single output frame.  

The particular format of the output files is specified by the user in the input config.ini file.  There are three options:
1)  full output - output all state variables; very large file sizes will results, but this is required for building initial conditions and for some analysis that require detailed composition and temperature information.  
2)  average state parameter output - species averaged temperature and velocity; electron density.  Probably best for most uses
3)  density only output - only electron density output.  Best for high-res instability runs where only the density is needed and the output cadence is high

The organization of the data in the MATLAB/octave workspace, after a single frame is loaded (via 'loadframe.m'), is as follows (MKSA units throughout):

### Time variables:  

simdate - a six element vector containing year, month, day, UT hour, UT minute, and UT seconds of the present frame

### Grid variables:  

<!--x1,x2,x3 - x1 is altitude (z in plots), x2 is east (x in plots), x3 north (y in plots); the sizes of these variables are stored in lxs by the MATLAB script.-->

structure xg - members xg.x1,2,3 are the position variables, xg.h\* are the metric factors, xg.dx\* are the finite differences, 

xg.glat,glon are the latitudes and longitudes (degrees geographic) of each grid point, xg.alt is the altitude of each grid point.  

xg.r,theta,phi - for each grid point:  radial distance (from ctr of Earth), magnetic colatitude (rads.), and magnetic longitude (rads.)

The grid structure, by itself, can be read in by the MATLAB function 'readgrid.m'; this is automatically invoked with 'loadframe.m' so there is not need to separately load the grid and output frame data.

### Temperature variable:

Ts (first three dimensions have size lxs; 4th dimension is species index:  1=O+,2=NO+,3=N2+,4=O2+,5=N+, 6=H+,7=e-)

### Density variable:  

ns (same indexing as temperature)

### Drifts:   

vs1 (same indexing as temperature)

x2-drift component:  v2 (same for all species, so this is just size lxs and is a 3D array)
x3-drift component:  v3

### Electromagnetic variables:  

current density:  J1, J2, J3
potential:  Phitop (EFL potential)

Note that the electric field is not included in the output file, but that it can be calculated from this output by taking -vxB at an altitude above about 200 km or by differentiating the top boundary electric potential 'Phitop' with respect to the x2 and x3 variables; however, note that if a curvilinear grid is used the derivatives must include the appropriate metric factors.


## Computing total electron content (TEC) and magnetic field perturbations

TEC and magnetic field variations can be calculated as a post-processing step in which the simulation data are read in and interpolated onto a regular geographic grid and then integrated accordingly using scripts in the './vis' directory - see 'TECcalc.m' and 'par\_magcalc.m'.  Note that script for computing magnetic fields, uses the MATLAB parallel processing toolbox, if available, to accelerate the calculations (which can take quite a while).  There is also a fortran parallel program for computing magnetic fields from large grids:  `magcalc.f90`.  An example of how to set this up is included in the `tohoku20113D_highres_var` examples in `./initialize`.


## References

The GEMINI3D model has been described and used in the following publications:  

Zettergren, M., & Semeter, J. (2012). Ionospheric plasma transport and loss in auroral downward current regions. Journal of Geophysical Research: Space Physics, 117(A6).

Zettergren, M. D., & Snively, J. B. (2013). Ionospheric signatures of acoustic waves generated by transient tropospheric forcing. Geophysical Research Letters, 40(20), 5345-5349.

Zettergren, M. D., & Snively, J. B. (2015). Ionospheric response to infrasonic‐acoustic waves generated by natural hazard events. Journal of Geophysical Research: Space Physics, 120(9), 8002-8024.

Zettergren, M. D., Semeter, J. L., & Dahlgren, H. (2015). Dynamics of density cavities generated by frictional heating: Formation, distortion, and instability. Geophysical Research Letters, 42(23).

Lynch, K. A., Hampton, D. L., Zettergren, M., Bekkeng, T. A., Conde, M., Fernandes, P. A., ... & Moen, J. (2015). MICA sounding rocket observations of conductivity‐gradient‐generated auroral ionospheric responses: Small‐scale structure with large‐scale drivers. Journal of Geophysical Research: Space Physics, 120(11), 9661-9682.

Perry, G. W., Dahlgren, H., Nicolls, M. J., Zettergren, M., St‐Maurice, J. P., Semeter, J. L., ... & Chen, S. (2015). Spatiotemporally resolved electrodynamic properties of a Sun‐aligned arc over Resolute Bay. Journal of Geophysical Research: Space Physics, 120(11), 9977-9987.

Fernandes, P. A., Lynch, K. A., Zettergren, M., Hampton, D. L., Bekkeng, T. A., Cohen, I. J., ... & Miceli, R. J. (2016). Measuring the seeds of ion outflow: Auroral sounding rocket observations of low‐altitude ion heating and circulation. Journal of Geophysical Research: Space Physics, 121(2), 1587-1607.

Swoboda, J., Semeter, J., Zettergren, M., & Erickson, P. J. (2017). Observability of ionospheric space‐time structure with ISR: A simulation study. Radio Science, 52(2), 215-234.

Zettergren, M. D., Snively, J. B., Komjathy, A., & Verkhoglyadova, O. P. (2017). Nonlinear ionospheric responses to large‐amplitude infrasonic‐acoustic waves generated by undersea earthquakes. Journal of Geophysical Research: Space Physics, 122(2), 2272-2291.

# The `magcalc.bin` Program

GEMINI includes a auxiliary program - `magcalc.f90` - that will compute magnetic field fluctuations (deviations from the Earth's main field) due to currents internal to by the model.  Magcalc reads in the output from a *completed* GEMINI disturbance simulation, namely the three components of current density and then uses the Biot-Savart (Ampere's) Law to compute magnetic fields from these currents by executing the appropriate volume integrals.

Magcalc will work for either 2D or 3D simulations; however, almost all testing has been conducted on 3D simulation *and* 2D simulations are likely inaccurate (will overestimate magnetic fluctuations usually) for purposes of computing field fluctuations from localized sources (longitudinally extended sources may still be okay).

Magcalc is fully parallelized, in particular making use of built-in mpi\_reduce functionality.  The full simulation grid is distributed to worker processes (domain parallelization) and each worker computes a piece of the Biot-Savart integral corresponding to their subdomain *for all field points* specified by user input.  The root process then collects integral portions from each worker and adds them together to form the full source grid Biot-Savart integral.

It is recommended that you run magcalc with the same number of processors and mpi image configuration as you ran the main simulation with.  This is not required but it will reduce the considerable computation time of the program.  Depending on the number of field points chosen by the user `magcalc.bin` can take as long or longer than the main program, particularly for a 3D set of field points.  For this reason, it is recommended that `magcalc.bin` be run with 1D or 2D lists of field points (e.g. ground plane or a line corresponding to spacecraft orbit) in order to reduce computation time.

It is worth noting that the use of the Biot-Savart law is subject to limitation, most notably we are not able to account of induced surface currents using this program.

## Running magcalc.bin

`magcalc.bin` can be invoked from the shell by:

```sh
mpirun -np 4 ./magcalc.bin /tmp/3d /tmp/3d/inputs/magfieldpoints.{h5,dat} -manual_grid <lid2> <lid3> <-debug>
```

This will compute magnetic fields over a grid at ground level using currents computed from a simulation stored in the directory `/tmp/3d`.  In order to run this program, you will need to create a set of field points at which the magnetic perturbations will be calculated - these are stored in the input file `magfieldpoints.{h5,dat}` used in the command above.  These could be a list of ground stations (irregular mesh), a regular mesh, or a set of satellite tracks (irregular mesh).  Optional inputs `lid2 lid3` are the number of mpi images to be used in the x2 and x3 directions, respectively, while using the `-debug` flag will cause the program to print a large amount of debug information to the console - this is useful for troubleshooting potential problems.  Note that the second argument listing the field point input file is optional so long as that file is named `magfieldpoints.{dat,h5}` and as long as that file is placed in `inputs/` subdirectory of the simulation directory.

There are additional command line argument parameters to magcalc that can be used to specify start and end times for magnetic field calculations; calculations are steady state so each time step is independent of all others (for a specified current density).  These options can invoked via:

```sh
mpirun -np 4 ./magcalc.bin /tmp/3d -start_time <YYYY> <MM> <DD> <UT_seconds> -end_time <YYYY> <MM> <DD> <UT_seconds>
```

The start and end time *must* correspond to output frames for the simulation.  If no start time is specified `magcalc.bin` will use the start time of the simulation; similarly, if no end time is specified `magcalc.bin` will use the simulation end time.  The files procressed by `magcalc.bin` will exclude the start time output file but include the end time output file.

## Creating Input Field Points for `magcalc.bin`

An example showing how to set up different types of input field point files for magcalc is shown here in the Moore, OK tornado example [https://github.com/gemini3d/gemci/tree/main/cfg/mag/arcs](https://github.com/gemini3d/gemci/tree/main/cfg/mag/arcs).

The input files for magcalc are organized as follows.  If raw binary input (.dat) is used then the input file contains:

```pseudo
1.  integer(4) :: lpoints                      ! the total number of field points at which the magnetic field is computed.
2.  real(8), dimension(lpoints) :: r,theta,phi ! arrays of spherical magnetic coordinates at which the magnetic field is to be computed
```

The above data must be present in addition to a variable `integer(4), dimension(3) :: gridsize` which indicates whether the list of points in the input file form a grid (this is useful for plotting routines which need to reshape the list/array into a proper multidimensional grid array.  If the input field points form a grid, the elements of `gridsize` are the number of grid points in the r,theta, and phi directions.  Otherwise the first element of gridsize is just lpoints, while the other two are -1, which indicates that the input points do *not* form a grid and should just be interpreted as a flat list of locations.


## Simluation vs. field point resolution

The Biot-Savart law involves both source locations (i.e. the grid the simulated currents are computed on) and field points (independent locations where we wish to evaluate the magnetic field) - see the mathematical formulation document description of the magnetic calculations for equations.  The source location resolution is given by the resolution at which the simulation of the currents has been run.

The field point grid can typically be *much* coarser particularly if you intend to evaluate currents on the ground - viz. away from the ionospheric source region.  This is a consequence of the fact that higher-order multipole moments (having smaller spatial structures) fall off quickly with distance from the current source such that they do not need to be resolved at far-field points.  Often a grid covering a ground range of mlats and mlons can be just 40 x 40 points and sufficiently resolve most structure.  If you intend to evaluate the magnetic fields in the ionosphere there will be quite a lot of small scale structure and you will want to use a larger grid, e.g. 192 x 192 or perhaps event large depending on the scale of the currents you wish to resolve.


## Output files created by `magcalc`

`magcalc` uses hdf5 output files containing the following variables:

```pseudo
real(wp), dimension(lpoints) :: Br, Btheta, Bphi.     ! there components of the magnetic field in up,south,east coordinates
```

Output files created by magcalc can be read using the mat_gemini interfact `gemini3d.read.magdata`.

## Visualizing Magnetic Field Perturbations Computed by `magcalc.bin`

The example script `magplot_fort_map.m` shows an example of how to load the results of running magcalc from a binary file into a MATLAB workspace and then plot these on a mapped grid.

<!--
One problematic aspect of magcalc is that you have to input the grid size into both the creation and plotting script and they must be consistent.  The corresponds to setting `ltheta` and `lphi` number of grid point in magnetic longitude and latitude in `gemini3d.model.magcalc()` and `gemini3d.plot.mag_map()`.  If these variables are not set propoerly the plotting program will not be able to read in, sort, and plot the data.  In the future this can be fixed by having magplot read in the grid size information from the input file that was created for the fortran program.
-->


## Caveats

`magcalc.bin` is not an entirely straightforward program to use and there are several issues that one must be aware of when running this program and analyzing the results.

### Cell Aspect Ratio

Due to the way the numerator and denominator in the Biot-Savart Law are averaged during the integration process some issues can arise when the resolutions in the different directions are vastly different.  Because distances are averaged over cells (which may be large in one dimension vs. another) there is a minimum distance which can be accomodated on a given grid; as such the field magnitudes can be underestimated in these cases.

The workaround for this issue is to avoid grids with highly anisotropic grid spacing; if this is not possible one could interpolate the model output into a more isotropic grid but there are no scripts currently that support this so the user would need to make these.

### Field Point Spacing

If one chooses a field point spacing *closer* than the source grid points, then some aliasing (spurrious periodic features) can occur.  Avoid such configurations as there is likely no advantage to having a more dense field point than source point grid.


### Parallel Currents are (likely) Required

When using magcalc it is best to have configured the GEMINI disturbance simulation to compute parallel currents by including the lines:

```nml
&Jpar
flagJpar=.true.
/
```

in the input configuration .nml file.  If you are using .ini file input then this parameter will be true by default and no changes are required.  While the program will work without a parallel current, it may produce results of dubious quality and `magcalc.bin` will list a warning in the console output.


<!--- MZ may add this later
## Example HPC queue submission script
-->

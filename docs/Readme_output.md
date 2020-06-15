# Gemini load and plot data

The default file format for Gemini is HDF5.
NetCDF4 file IO is [optionally available](./Readme_cmake.md).
Raw file output was the original Gemini format, but does not support all features, is not tested and not recommended.


## Default plotting

To produce the default style plots for simulation run (python):

```sh
gemini_plot /tmp/mysim
```

The MATLAB version of this script can be invoked by starting MATLAB and running:

```matlab
direc='~/simulations/mysim'
plotall(direc,{'png','eps'})
```
where the latter argument specifies the file type to which matlab will print the plots.

The script reads a sequence of files from a simulation.
This script saves a copy of the output plots into the simulation output directory.


## Loading simulation data

To load output data from a GEMINI simulation, e.g. for purposes of custom plotting or additional calculation, a ```loadframe``` API is provided.  Loadframe APIs work with any type of GEMINI file output, viz. with binary, hdf5, or netcdf output.

In python:

```sh
MH - please add commands to load variables into python workspace or whatever it is called.
```

In MATLAB:

```matlab
dat=loadframe(direc,ymd,UTsec,flagoutput,mloc,xg,config_file,realbits)
```
This will load into the structure dat all of the plasma information from the simulation frame corresponding to output directory ```direc``` on date (year, month, day) ```ymd``` at time (UT seconds) ```UTsec```.  The remaining arguments are optional and are mainly present to prevent the code from having to reload, e.g. the grid, if it is called repeatedly for different time frames in the same simulation.  The output dat contains fields ```dat.ne,dat.Te,dat.Ti,dat.v1``` etc. corresponding to different calculated variables of interest (see descriptions below).

By default the ```loadframe``` API will not load the grid or simulation configuration information.

To load the simulation information:

In Python:

```sh
MH - please add commands to load variables into python workspace or whatever it is called.
```

In MATLAB:

```matlab
[ymd,UTsec,tdur,dtout,flagoutput,mloc,activ,indat_size,indat_grid,indat_file] = readconfig(path)
```

To load the grid a ```readgrid``` function is provided:

In Python:

```python
MH - please add commands to load variables into python workspace or whatever it is called.
```

In MATLAB

```matlab
xg=readgrid(path,realbits);
```
Here ```path``` is the path to the grid data file and ```realbits``` indicates 32-bit or 64-bit floating point data in the file.


## Output file content

The particular format of the output files is specified by the user in the input config.nml file.  There are three options:

1. full output (```flagoutput=1``` in the config.nml file) - output all state variables; very large file sizes will results, but this is required for building initial conditions and for some analysis that require detailed composition and temperature information.  Selecting full output will require a large amount of disk space  - *up to several terabytes per simulation*.  The other use for full output is that it is needed for any kind of simulation milestone/restart functionality.
2. average state parameter output (```flagoutput=2```) - species averaged temperature and velocity; electron density.  Probably best for most uses not requiring a full output of the plasma state.
3. density only output (```flagoutput=3```) - only electron density output.  Best for high-res instability runs where only the density is needed and the output cadence is high.  This choice for output saves an enormous amount of disk space.

MKSA units are used throughout.


## Time variables

```simdate``` - a six element vector containing year, month, day, UT hour, UT minute, and UT seconds of the present frame


## Grid variables

<!--x1,x2,x3 - x1 is altitude (z in plots), x2 is east (x in plots), x3 north (y in plots); the sizes of these variables are stored in lxs.-->

All three dimensional arrays have dimensions ordered as ```x1,x2,x3```, i.e. the first dimension corresponds to the *field-line* coordinate.

structure xg - members ```xg.x1,2,3``` are the position variables, ```xg.h*``` are the metric factors, ```xg.dx*``` are the finite differences.  For Cartesian coordinates x1 is altitude, x2 is eastward distance, x3 is northward distance (all meters).

```xg.glat,glon``` are the latitudes and longitudes (degrees geographic) of each grid point, ```xg.alt``` is the altitude of each grid point.

```xg.r,theta,phi``` - for each grid point:  radial distance (from ctr of Earth), magnetic colatitude (rads.), and magnetic longitude (rads.).  The magnetic pole and moment is hard coded into the grid generation scripts.

### Temperature variable

```Ts``` (first three dimensions have size lxs; 4th dimension is species index:  1=O+,2=NO+,3=N2+,4=O2+,5=N+, 6=H+,7=e-)

### Density variable

```ns``` (same indexing as temperature)

### Parallel to **B** (x1) Drifts

```vs1``` (same indexing as temperature)

x2-drift component:  ```v2``` (same for all species, so this is just size lxs and is a 3D array)
x3-drift component:  ```v3```

### Electromagnetic variables

current density:  ```J1, J2, J3```
potential:  ```Phitop``` (EFL potential)

Note that the electric field is not included in the output file, but that it can be calculated from this output by taking -vxB at an altitude above about 200 km or by differentiating the top boundary electric potential 'Phitop' with respect to the x2 and x3 variables; however, note that if a curvilinear grid is used the derivatives must include the appropriate metric factors.


## Computing total electron content (TEC)

TEC and magnetic field variations can be calculated as a post-processing step in which the simulation data are read in and interpolated onto a regular geographic grid and then integrated accordingly using scripts in the './vis' directory - see `TECcalc.m`.
An example of how to plot TEC computed by this script is included in `TECplot_map.m` in Gemini-matlab repo


## Computing magnetic field perturbations

Magnetic field perturbations for a given simulation can be computed, after the fact, using a separate fortran mpi program, magcalc.  Detailed instructions are included in the magcalc readme, [./Readme_magcalc.md](./Readme_magcalc.md).

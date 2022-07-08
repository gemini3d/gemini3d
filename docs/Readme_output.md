# Gemini load and plot data

The file format for Gemini is HDF5.

## Default plotting

To produce the default style plots for simulation run (python):

```python
import gemini3d.plot
direc="~/simulations/mysim"
gemini3d.plot.plot_all(direc,saveplot_fmt="png")
```
The last argument is the format in which the plots will be saved, either "png" or "eps".

The MATLAB version of this script can be invoked by starting MATLAB and running:

```matlab
run <location of matgemini>/setup.m
direc='~/simulations/mysim'
gemini3d.plot.plotall(direc,"png")
```
where the latter argument specifies the file type to which matlab will print the plots.  Either "png" or "eps" are currently supported.

The Python/MATLAB ```plotall``` script reads a sequence of files from the simulation directory corresponding to the full sequence of outputs.  This script saves a copy of the output plots into the simulation output directory under "./plots/".


## Loading simulation data

To load output data from a GEMINI simulation, e.g. for purposes of custom plotting or additional calculation, a `loadframe` API is provided.  Loadframe APIs work with any type of GEMINI file output, viz. with HDF5 output.

In python:

```sh
import gemini3d.read
direc="~/simulations/mysim"
dat=gemini3d.read.frame(direc,time=<datetime_variable>)
```

If you have already read in the simulation config info (see below) often this will take the form:

```sh
import gemini3d.read
direc="~/simulations/mysim"
cfg = gemini3d.read.config(direc)
dat=gemini3d.read.frame(direc,time=cfg["time"][-1])
```

The output dictionary ```dat``` will have xarray entries corresponding to simulation variables corresponding to the output frame requested.

In MATLAB:

```matlab
time=datetime([year,month,day,hour,minute,second])
dat = gemini3d.loadframe(direc, "time", time)
```
This will load into the structure dat all of the plasma information from the simulation frame corresponding to output directory `direc` on time=datetime().  The time variable must be a MATLAB datetime structure as shown above.  The remaining arguments are optional and are mainly present to prevent the code from having to reloadn data if it is called repeatedly for different time frames in the same simulation.  The output dat contains fields `dat.ne,dat.Te,dat.Ti,dat.v1` etc. corresponding to different calculated variables of interest (see descriptions below).


## Loading simulation metadata

By default the `loadframe` API will not load the grid or simulation configuration information; these require separate calls:  To load the simulation configuration information in Python:

```Python
import gemini3d.read
direc="~/simulations/mysim"
cfg=gemini3d.read.config(direc)
```

The ```cfg``` dictionary will contain all of the information from the config.nml or config.ini file.

In MATLAB the config file can be read via:

```matlab
run <location of matgemini>/setup.m
direc='~/simulations/mysim'
cfg = gemini3d.read.config(direc)
```

The returned ```cfg``` variable is a MATLAB struct and will contain input file information.


## Loading simulation grid

To load the grid a `readgrid` function is provided:

In Python:

```python
import gemini3d.read
direc="~/simulations/mysim"
xg=gemini3d.read.grid(direc)
```

In MATLAB

```matlab
run <location of matgemini>/setup.m
direc='~/simulations/mysim'
xg = gemini3d.read.grid(direc)
```
Here ```direc``` is the path to the grid data file.  And the output object will be a dictionary or structure containing all of the grid information, list Readme_input.md for a full list.


## Output file content

The particular format of the output files is specified by the user in the input config.nml file.  There are three options:

1. full output (`flagoutput=1` in the config.nml file) - output all state variables; very large file sizes will results, but this is required for building initial conditions and for some analysis that require detailed composition and temperature information.  Selecting full output will require a large amount of disk space  - *up to several terabytes per simulation*.  The other use for full output is that it is needed for any kind of simulation milestone/restart functionality.
2. average state parameter output (`flagoutput=2`) - species averaged temperature and velocity; electron density.  Probably best for most uses not requiring a full output of the plasma state.
3. density only output (`flagoutput=3`) - only electron density output.  Best for high-res instability runs where only the density is needed and the output cadence is high.  This choice for output saves an enormous amount of disk space.

MKSA units are used throughout.


## Time variables

`simdate` - a six element vector containing year, month, day, UT hour, UT minute, and UT seconds of the present frame


## Temperature variable

`Ts` (first three dimensions have size lxs; 4th dimension is species index:  1=O+,2=NO+,3=N2+,4=O2+,5=N+, 6=H+,7=e-)

## Density variable

`ns` (same indexing as temperature)

## Parallel to **B** (x1) Drifts

`vs1` (same indexing as temperature)

x2-drift component:  `v2` (same for all species, so this is just size lxs and is a 3D array)
x3-drift component:  `v3`

## Electromagnetic variables

current density:  `J1, J2, J3`
potential:  `Phitop` (EFL potential)

Note that the electric field is not included in the output file, but that it can be calculated from this output by taking -vxB at an altitude above about 200 km or by differentiating the top boundary electric potential 'Phitop' with respect to the x2 and x3 variables; however, note that if a curvilinear grid is used the derivatives must include the appropriate metric factors.


## Computing total electron content (TEC)

TEC and magnetic field variations can be calculated as a post-processing step in which the simulation data are read in and interpolated onto a regular geographic grid and then integrated accordingly using scripts in the './vis' directory - see `TECcalc.m`.
An example of how to plot TEC computed by this script is included in `TECplot_map.m` in Gemini-matlab repo


## Computing magnetic field perturbations

Magnetic field perturbations for a given simulation can be computed, after the fact, using a separate fortran mpi program, magcalc.
Detailed instructions are included in the magcalc readme, [./Readme_magcalc.md](./Readme_magcalc.md).

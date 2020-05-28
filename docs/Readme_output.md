# Gemini load and plot data

The default file format for Gemini is HDF5.
NetCDF4 file IO is [optionally available](./Readme_cmake.md).
Raw file output was the original Gemini format, but does not support all features, is not tested and not recommended.

## Plot

```sh
python plotall.py /tmp/mysim
```

The script reads a sequence of files from a simulation.
This script saves a copy of the output plots into the simulation output directory.

The particular format of the output files is specified by the user in the input config.nml file.  There are three options:

1. full output - output all state variables; very large file sizes will results, but this is required for building initial conditions and for some analysis that require detailed composition and temperature information.
2. average state parameter output - species averaged temperature and velocity; electron density.  Probably best for most uses
3. density only output - only electron density output.  Best for high-res instability runs where only the density is needed and the output cadence is high

MKSA units are used throughout.

## Time variables

simdate - a six element vector containing year, month, day, UT hour, UT minute, and UT seconds of the present frame

## Grid variables

<!--x1,x2,x3 - x1 is altitude (z in plots), x2 is east (x in plots), x3 north (y in plots); the sizes of these variables are stored in lxs.-->

structure xg - members xg.x1,2,3 are the position variables, `xg.h*` are the metric factors, `xg.dx*` are the finite differences,

xg.glat,glon are the latitudes and longitudes (degrees geographic) of each grid point, xg.alt is the altitude of each grid point.

xg.r,theta,phi - for each grid point:  radial distance (from ctr of Earth), magnetic colatitude (rads.), and magnetic longitude (rads.)

### Temperature variable

Ts (first three dimensions have size lxs; 4th dimension is species index:  1=O+,2=NO+,3=N2+,4=O2+,5=N+, 6=H+,7=e-)

### Density variable

ns (same indexing as temperature)

### Drifts

vs1 (same indexing as temperature)

x2-drift component:  v2 (same for all species, so this is just size lxs and is a 3D array)
x3-drift component:  v3

### Electromagnetic variables

current density:  J1, J2, J3
potential:  Phitop (EFL potential)

Note that the electric field is not included in the output file, but that it can be calculated from this output by taking -vxB at an altitude above about 200 km or by differentiating the top boundary electric potential 'Phitop' with respect to the x2 and x3 variables; however, note that if a curvilinear grid is used the derivatives must include the appropriate metric factors.

## Computing total electron content (TEC)

TEC and magnetic field variations can be calculated as a post-processing step in which the simulation data are read in and interpolated onto a regular geographic grid and then integrated accordingly using scripts in the './vis' directory - see `TECcalc.m`.
An example of how to plot TEC computed by this script is included in `TECplot_map.m` in Gemini-matlab repo


## Computing magnetic field perturbations

Magnetic field perturbations for a given simulation can be computed, after the fact, using a separate fortran mpi program, magcalc.  Detailed instructions are included in the magcalc readme, [./Readme_magcalc.md](./Readme_magcalc.md).
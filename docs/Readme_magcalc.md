# The magcalc.bin fortran program

GEMINI includes a auxiliary program - ```magcalc.f90``` - that will compute magnetic field fluctuations (deviations from the Earth's main field) due to currents computed by the model.  Magcalc reads in the output from a *completed* GEMINI disturbance simulation, namely the three components of current density and then uses the Biot-Savart (Ampere's) Law to compute magnetic fields directly from these currents.  Magcalc will work for either 2D or 3D simulations.   

Magcalc is fully parallelized, in particular making use of built-in mpi *reduce* functionality.  The full simulation grid is distributed to worker processes (domain parallelization) and each worker compute the piece of the Biot-Savart integral corresponding to their subdomains and all field points at which the field is to be evaluated (as specified by user input).  The root process then collects integral portions from each worker and adds them together to form an estimate of the Biot-Savart integral over the entire grid.  

It is recommended that you run magcalc with the same number of processors and mpi image configuration as you ran the main simulation with.  


## Caveats

When using magcalc it is important to insure you have configured the GEMINI disturbance simulation to compute parallel currents by including the lines: 

```
&Jpar
flagJpar=.true.
/
```
in the input configuration .nml file.  If you are using .ini file input then this parameter will be true by default and no changes are required.  While the program will work without a parallel current it may produce results of dubious quality.  


## Running magcalc.bin

Note that there is also a utility that can compute magnetic fields from the currents calculated by GEMINI.
This can be run by:

```
mpirun -np 4 ./magcalc.bin /tmp/3d tests/data/test3d/input/magfieldpoints.dat <lid2 lid3> <-debug>
```

This will compute magnetic fields over a grid at ground level using currents computed from a simulation stored in the directory ```/tmp/3d```.  In order to run this program, you will need to create a set of field points at which the magnetic perturbations will be calculated - these are stored in the input file ```magfieldpoints.dat``` used in the command above.  These could be a list of ground stations (irregular mesh), a regular mesh, or a set of satellite tracks (irregular mesh).  Optional inputs ```lid2 lid3``` are the number of mpi images to be used in the x2 and x3 directions, respectively, while using the ```-debug``` flag will cause the program to print a large amount of debug information to the console - this is useful for troubleshooting potential problems.   

Magcalc is a computationally intensive program and it is suggested that you run it with the same number of mpi images as the corresponding disturbance simulations.  If a large grid was used with your simulation, magcalc will consume a large amount of memory, as well (commensurate with memory use by the base GEMINI program).  


## Creating input field points for ```magcalc```

An example showing how to set up different types of input field point files for magcalc is shown here in the Moore, OK tornado example [https://github.com/gemini3d/gemini-examples/tree/master/initialize/mooreOK3D_hemis](https://github.com/gemini3d/gemini-examples/tree/master/initialize/mooreOK3D_hemis).  

Other examples using magcalc setup scripts include:

1. Iowa thunderstorm [https://github.com/gemini3d/gemini-examples/tree/master/initialize/iowa3D](https://github.com/gemini3d/gemini-examples/tree/master/initialize/iowa3D).
2. Tohoku earthquake (low resolution) [https://github.com/gemini3d/gemini-examples/tree/master/initialize/tohoku20113D_lowres](https://github.com/gemini3d/gemini-examples/tree/master/initialize/tohoku20113D_lowres).

The low-resolution Tohoku examples, in particular, is one of the cheapest simulations to run (will run on a laptop in several hours) that will provide a good context in which to compute magnetic field perturbations.  


## Visualizing magnetic field perturbations computed by magcalc.f90

The example script `magplot_fort_map.m` shows an example of how to load the results of running magcalc from a binary file into a MATLAB workspace and then plot these on a mapped grid.

One problematic aspect of magcalc is that you have to input the grid size into both the creation and plotting script and they must be consistent.  The corresponds to setting ```ltheta``` and ```lphi``` number of grid point in magnetic longitude and latitude in both the ```magcalc_setup.m``` scripts and ```magplot_fort_map.m``` scripts.  If these variables are not set propoerly the plotting program will not be able to read in, sort, and plot the data.  In the future this can be fixed by having magplot read in the grid size information from the input file that was created for the fortran program.  


<!--- MZ may add this later
## Example HPC queue submission script

```

```
-->

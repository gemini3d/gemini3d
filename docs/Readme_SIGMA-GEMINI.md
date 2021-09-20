# Setup and run SIGMA-GEMINI interface simulations - VEGA version

Before running the interface, make sure that the GEMINI is setup, built and running as mentioned in the link below

https://github.com/gemini3d/gemini3d/blob/main/docs/Readme_VEGA.md

The GEMINI outputs the simulation results for the GEMINI-example that user wants to run (here we use "arcs" as an example). Note the path of the output files.

The next step is to map the GEMINI output as an input to the SIGMA using a MATLAB interface function (SIGMAinterp_time_MZfunction)

Open the MATLAB file and edit the path of the directory to the location of GEMINI output files that user wants to map.

```PowerShell
Example: 
direc = ('~/Projects/GEMINI/arcs/');
```
#### Read the configuration and grid data

User needs to read the configurationa and grid data from the GEMINI output files. This helps to extract the GEMINI configuration data (e.g. GEMINI data output time) and the grid data (e.g. grid dimensions).

```PowerShell
Example:
cfg = gemini3d.read.config(direc)  % configuration data
xg = gemini3d.read.grid(direc)     % grid data
```
#### Extract the data

```PowerShell
Example:
dtout = cfg.dtout                  % differential time
dateval = cfg.times                % Start and end times
x1 = xg.x1                         % grid dimension in x1 direction (Similarly for x2 and x3 directions)
```

#### Load the simulation data

Here the plasma paramters are read and loaded into the function. User can choose to load either all the GEMINI plasma parameters or the required parameter. 
```PowerShell
Example:
datplasma = gemini3d.read.frame(direc,"time",dateval,"vars","ne")
```

Once the plasma data is read into the interface function, it interpolates the data in time and sends out the interpolated plasma data to the main SIGMA function.


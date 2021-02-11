# Test MAGIC-GEMINI runs

This section demonstrated the pathways to create various MAGIC (neutral gas model) inputs to GEMINI simulations.

1. [MAGIC simulation description](#magic_description)
2. [3D Dipole GEMINI simulation with 3D Cartesian MAGIC input](#3d_dipole_3d_neutral)
3. [2D Dipole GEMINI simulation with 2D Cartesian MAGIC input](#2d_dipole_2d_neutral)
4. [2D/3D Dipole GEMINI simulations with 2D Axisymmetric MAGIC input](#3d_dipole_2daxi_neutral)

<a name="magic_description"></a>
## 1. MAGIC simulation description

Neutral atmospheric model MAGIC simulation is run with next configuration:

1. Type of the domain: 3D Cartesian;
2. Physical domain size: 750 x 750 x 500 km (meridional x zonal x vertical);
3. Numerical domain size: 90 x 80 x 100 points (in meridional x zonal x vertical directions);
4. Resolution: 8333.3 x 9375.0 x 5000 m (in meridional x zonal x vertical directions);
5. Simulation time: 1800 sec with 5 sec time step for output.

Main specification of the simulation:

```
90              mx          = cells in x direction
80              my          = cells in y direction
100             mz          = cells in z direction
360             nout        = number of output times to print results
1               outstyle    = style of specifying output times
1800            tfinal      = final time
1               dtv(1)	  = initial dt (used in all steps if method(1)=0)
0.0d0           t0          = initial time
-3.75d5         xlower	 = left edge of computational domain
3.75d5          xupper      = right edge of computational domain
-3.75d5         ylower	 = bottom edge of computational domain
3.75d5          yupper	= top edge of computational domain
0.0d0           zlower      = front edge of computational domain
5.0d5           zupper      = back edge of computational domain
```

The source function represents a slightly moving (in zonal direction from the center of the domain) axisymmetric source:

```
 1        forcemth     = 0 for no source, 1 for one, 2 for 'source.data' file
 0.035    omega	= forcing frequency
 1.5	   amplitude    = forcing magnitude (m/s or kg*m/s^2)
 0	      propx        = horizontal-x propagation constant
 2d-4	   propy        = horizontal-y propagation constant
 0        xpos         = x-axis position of oscillator
 0        ypos         = y-axis position of oscillator
 0        zpos         = z-axis position of oscillator
 2d4	   xwidth	= x-streamwise "width" of gaussian envelope
 2d4	   ywidth	= y-spanwise "width" of gaussian envelope
 2e3	   zwidth	= x-vertical "height" of gaussian envelope
 180	   tcenter      = peak forcing time (center of gaussian packet)
 180	   twidth	= temporal "width" of gaussian envelope
 0        vsrcx        = x-axis oscillator velocity
```

Simulation outputs are stored in HDF5 files separately for each time step (the data is stored [here](https://www.dropbox.com/sh/z448k1ttjm9ljvx/AADKSWqPtiMuKYxYDmrw8kTHa?dl=0). The animation of fluid dynamics in 3 direction is [here](https://www.dropbox.com/s/pt4lqavc717a90a/MAGIC-movie.mp4?dl=0) and can be replicated using [animation.m](https://www.dropbox.com/s/pt4lqavc717a90a/MAGIC-movie.mp4?dl=0) from the MAGIC simulation folder:

![Screenshot from MAGIC simulation](pics/magicexample.jpg)

<a name="3d_dipole_3d_neutral"></a>
## 2. 3D Dipole GEMINI simulation with 3D Cartesian MAGIC input

GEMINI grid is 3D Dipole of the sizeTo prepare GEMINI neutral inputs, run Matlab script ``initialize3DCARD.m`` (``readandoutput3DCART.m`` required).
GEMINI neutral input requires the specification of simulation initial time (UT), time step (sec) and grid structure, as shown below:

```
% Initial time of GEMINI simulation and time step of neutral inputs
ymd0=[2011,3,11];
UTsec0=35100.0+5;
dtneu=5;

% Specify input data
lt=361; % number of time steps
lx1=80; % number of points in zonal direction
lx2=90; % number of points in meridional direction
lx3=100; % number of points in vertical direction
```

Note, that ``lx1`` is zonal, ``lx2`` - meridional, and ``lx3`` - vertical directions, whereas the structure of matrix to be saved as GEMINI input is ``[vertical,zonal,meridional]``. Matrix permutation should be done accordingly in case of need. The figure below is a snapshot from the GEMINI simulation with specified 3D Cartesian MAGIC source. Simulation data and figures can be found [here](https://www.dropbox.com/sh/59s2v9boscs61bb/AADseeIShD2tiFHlnaeIgFcNa?dl=0).

![GEMINI simulation output with 3D Cartesian MAGIC input](pics/3Dcart-v1-20110311_36225.000000.png)

<a name="2d_dipole_2d_neutral"></a>
## 3. 2D Dipole GEMINI simulation with 2D Cartesian MAGIC input

To prepare GEMINI neutral inputs, run Matlab script ``initialize3DCARD.m`` (``readandoutput3DCART.m`` required).
GEMINI neutral input requires the specification of simulation initial time (UT), time step (sec) and grid structure, as shown below:

```
% Initial time of GEMINI simulation and time step of neutral inputs
ymd0=[2011,3,11];
UTsec0=35100.0+5*Frame;
dtneu=5;

% Specify input data
lt=361; % number of time steps
lx1=90; % latitude
lx2=100; % altitude
```

Note, that here ``lx1`` is meridional, and ``lx2`` - vertical directions, whereas the structure of matrix to be saved as GEMINI neutral input is ``[vertical,meridional]``. Matrix permutation should be done accordingly in case of need. The figure below is a snapshot from the GEMINI simulation with specified 2D Cartesian MAGIC source. Simulation data and figures can be found [here](https://www.dropbox.com/sh/74gvxzoqyf459do/AAC92yoEkDH0CHXhX5CPmV4va?dl=0).

![GEMINI simulation output with 3D Cartesian MAGIC input](pics/2Dcart-v1-20110311_36225.000000.png)

<a name="3d_dipole_2daxi_neutral"></a>
## 4. 2D/3D Dipole GEMINI simulation with 2D Axisymmetric MAGIC input

To prepare GEMINI neutral inputs, run Matlab script ``initialize3DCARD.m`` (``readandoutput3DCART.m`` required).
GEMINI neutral input requires the specification of simulation initial time (UT), time step (sec) and grid structure, as shown below:

```
% Initial time of GEMINI simulation and time step of neutral inputs
ymd0=[2011,3,11];
UTsec0=35100.0+5*Frame;
dtneu=5;

% Specify input data
lt=361; % number of time steps
lx1=45; % radial direction
lx2=100; % upward direction
```

Note, that here ``lx1`` is a radial, and ``lx2`` - vertical directions, whereas the structure of matrix to be saved as GEMINI neutral input is ``[vertical,radial]``. Matrix permutation should be done accordingly in case of need. Be sure that point 1 of matric is radial direction corresponds to the r=0. For reference, below is a snapshot from MAGIC matrix from Matlab ready to output into GEMINI neutral file of vertical fluid velocity ``pcolor(velzfull)`` (note that in depreciated *.dat MAGIC inputs, the matrix should be reversed):

![Reference for Axisymmetric matrix structure](pics/pcolorreference.jpg)

The same MAGIC inputs can be used for both 2D and 3D GEMINI simulations. The figures below are a snapshot from the 2D and 3D GEMINI simulation with specified 2D Axisymmetric MAGIC source. Simulation data and figures can be found [here](https://www.dropbox.com/sh/m5juo70c9n3j0ud/AAC-lNyltlBL7sbzWOq0L0A6a?dl=0) (3D) and [here](https://www.dropbox.com/sh/oi7f3i5lf2bxtbr/AAAEpT6MegGXxeQ3MRUyvnUFa?dl=0) (2D).

2D Dipole GEMINI simulation with Axisymmetric MAGIC source:
![GEMINI simulation output with 3D Cartesian MAGIC input](pics/2daxi-v1-20110311_36225.000000.png)

3D Dipole GEMINI simulation with Axisymmetric MAGIC source:
![GEMINI simulation output with 3D Cartesian MAGIC input](pics/3daxi-v1-20110311_36225.000000.png)

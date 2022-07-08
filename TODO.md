# Development Roadmap and "To Do" List

This file is intended to document development priorities for the GEMINI project.


## Curated Examples to Develop

* Equatorial plasma bubble simulation examples and tweaks
* Simulation of Perkins instability (this was started but needs to be merged with master)
* Global scale (mid/low) modeling example, need to include EIA per Paul Inchin, need to flag periodic in x3 variable (longitude)


## Changes to Mathematical/Physical Formulation of GEMINI

### High Priority

* Heat flux boundary conditions for the electrons (possibly useful for SAID/STEVE)

### Requiring Resources

* Small-scale simulation development for simulation of turbulence:
	* Diamagnetic drift and perpendicular ambipolar fields - necessary for the smallest scales, e.g. less than 100 m
	* Resolved potential solutions - decimate parallel grid down to Farley mapping scale for perp resolution then so the solve on that coarse grid then interpolate back up to original grid.  I've had luck with MUMPS solves in reasonable time up to 300 x 300 x 15 grid points which is probably enough to do something interesting with appropriate periodic and Lagrangian grids (moving at E x B).
	*  For smallest scale simulations it makes sense to have the simulation able to use a Lagrangian frame of reference to reduce the total number of grid points needed.
* Option for true coordinates to be used in the computations of magnetic perturbations (instead of flattened-out spherical)
* Adaptive Mesh Refinement
	* p4est based; ForestClaw?
* Two way coupling with MAGIC
	* High-latitude is the obvious application for this
* GLOW application on dipole geomagnetic field lines
* Equatorial spread F simulation, e.g. need to add gravitation drift terms and source terms.


## Refactoring and/or Cleanup Needed

### High Priority

* Magnetic field calculation code:
	* Mag field points need a separate file size so that it isn't hard coded in scripts...  This should be handled similarly to the neutral inputs and precipitaiton/potential boundary inputs...
	* Redo the way magcalc and magplot deal with file names and locations also have the shape of the grid determined automatically via some sort of input file.
	* magplot needs some way to know what the grid dimensions are so the user doesn't have to manually input
* Fair bit of code repetition in top-level electric field and precipitation interpolation routines
* Boundary condition modules for the electrodynamics and precipitation should be removed in favor of submodules of the electrodynamics and ionization modules.  If we do this are we breaking backwards compatibility with older compilers?  Do we even care?
* MSISmatlab is a mess, uses dmy instead of ymd and UThrs instead of UTsec - may break some other scripts if we change this.

### Lower Priority

* Remove the array permuting form the fortran code and do this from the MATLAB/octave scripts.  These scripts should provide a permutation array for the dimensions to the fortran code (e.g. [1,2,3] or [3,1,2], which are even, or [1,3,2] which is odd), which then knows if the coordinate system is right-handed or left-handed so it can adjust the cross products accordingly.

### Completed

* (Completed) Code duplication in electrodynamics module (haloing part should be written as a subroutine).  This duplication also exists in other places, e.g., in RK2prep - this has now become a serious readability issue for anyone trying to modify those files
* (Completed) Axisymmetric and Cartesian interpolations should be combined (much code-sharing)
* (Completed) There are now numerous versions of routines corresponding to message passing in x3 vs. on a x2/x3 process grid.  Somehow the x3 routines need to be kept as they may be faster in some (hopefully unusual) situations.  Michael suggests a submodule...
* (Completed) Some modules have now become excessively large, e.g. mpimod and calculus...  These need to be organized and split up
* Handling of metric factors in the potential solves is sloppy - need to be passing into solver and used to eval. geometric terms there - would be more clear to reader...
* Elliptic solvers do not need to check for root vs. workers anymore; is done from calling functions
* (Completed) Clean up neutral interpolation code


## Simulation Milestones and Restarting

* Fix restart code for precipitation and electric field input files (need to be primed like the neutral input)
* Proper restarting will require reading in an initial potential value, as well - means that the matlab input scripts need to be fixed/updated too...
* Milestone full-data output to enable restarting, need some metadata in output in order for code to know when the milestone outputs are done


## Grid Generation and I/O

* Magnetic pole location maybe should be taken as input to geo*2geo*2, needs to be updated or to allow user-defined secular variations depending on epoch, etc.
* Magnetic moment should not be hardcoded, inconsistent with IGRF, maybe need a lookup table based on latest IGRF.
* Chunk input data directly from file to workers.


## Issues affecting ease of use

* Work on some means for plot color axes to be adjustable.
* Refactor the precip and Efield generation scripts (possibly already done by Michael H.).  Document exactly the expectations that the fortran code has for these input files.
* (low priority) Script that will run a simulation and then plot the results under the same session


## Continuous Integration

### High Priority

* Unit tests for hyperbolic solver

### Completed

* Unit test for parabolic solver
* Unit test for basic elliptic solver


## Coding style and standards issues

* There may be a performance boost by using the Fortran 2008 `contiguous` attribute on the `pointer` arrays where right now it is manually repacked--`contiguous` means we DON'T repack manually, the compiler will repack IF and ONLY IF it needs too.  We may get a performance boost by eliminating manual repacking and using `real, contiguous, pointer` instead. [Reference](https://modelingguru.nasa.gov/servlet/JiveServlet/previewBody/1527-102-1-2631/N1729-4.pdf) page 7.


## Feature requests

### High Priority

* Ability to specify a background neutral atmospheric wind profile that is constant over the entire grid.
* Include precipitation information (Qp,E0p,mlatp,mlonp variables from precipBCs module) in GEMINI output files?
* (EFFICIENCY) Exclusion of null points from field aligned advection, thermal conduction, and source terms - could improve performance
* (Underway) Option to run the code in a single precision mode - would help with memory limited systems although it's not clear how this would impact numerics (I've never tested my methods in single precision)

### Completed

* (Completed) Parallel domain decomposition in x2 *and* x3 - this is a big task that is likely to be left aside until I can renew funding.  It's also questionable how useful it is at this point where my typical runs are 32-256 cores (although undoubtedly it may become useful for runs with thousands of cores).  I've found good speedup even dividing the x3 dimension into slabs 2 grid points wide; although that means passing essentially all the grid data around via mpi, the large number of operations per slab means that the effective overhead here is not too much to prevent this from being useful.  However, for simulations that run with GLOW this will massively speed things up...
* (Completed) Add 3Dtest to ctest
* (Completed) HDF5 file input and output
* (Completed) Periodically updating background neutral atmosphere - should really be done for simulations more than a few hours long but will affect performance

## Interfaces with other models

### High Priority

* Pass inclination from GEMINI to GLOW rather than using IGRF in glow while using a dipole in GEMINI
* Return superthermal current from GLOW in order to calculate the thermal current in GEMINI
* Inverted grid must be passed to GLOW if running a curvilinear altitude array to GLOW
* Talk to GLOW developers about how GLOW might be used on closed field lines for GEMINI (not sure if it can be used in this way currently)
* Have GLOW output VER in order to deal with weird observing geometries.

### Completed

* (Completed) Ability to use GLOW to compute ionization and heating rates, as well as brightnesses of various bands of interest.  This brings up a lot of questions about how GLOW will function on a closed field-line grid; we may need to talk to Stan about this.

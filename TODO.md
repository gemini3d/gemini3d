# Development roadmap and to do list

This file is intended to document development priorities for the GEMINI project.  


## Parallel in x2 and x3 changes
* grid.f90 needs to create x2all
* derivatives must know about x2 division
* haloing of divs must know about x2 division
* rewrite custom broadcast and gather ops
<!--* alias halo functions -->


## Things that need to be checked or otherwise dealt with

* Work on some means for plot color axes to be adjustable.
* refactor the precip and Efield generation scripts
* Redo the way magcalc and magplot deal with file names and locations
* magplot needs some way to know what the grid dimensions are so the user doesn't have to manually input 
* Fix restart code for precipitation and electric field input files (need to be primed like the neutral input)
* run + plot script...


## Future Code refactoring

* Fair bit of code repetition in top-level electric field and precipitation interpolation routines
* Code duplication in electrodynamics module (haloing part should be written as a subroutine)
* Cleanup of BCs interpolation source files
* Axisymmetric and Cartesian interpolations should be combined (much code-sharing)
* Remove the array permuting form the fortran code and do this from the MATLAB/octave scripts.  These scripts should provide a permutation array for the dimensions to the fortran code (e.g. [1,2,3] or [3,1,2], which are even, or [1,3,2] which is odd), which then knows if the coordinate system is right-handed or left-handed so it can adjust the cross products accordingly.
* Boundary condition modules for the electrodynamics and precipitation should be removed in favor of submodules of the electrodynamics and ionization modules.  If we do this are we breaking backwards compatibility with older compilers?  Do we even care?
* MSISmatlab is a mess, uses dmy instead of ymd and UThrs instead of UTsec - not sure what this will affect is we change...


## Coding style and standards issues
* There may be a performance boost by using the Fortran 2008 `contiguous` attribute on the `pointer` arrays where right now it is manually repacked--`contiguous` means we DON'T repack manually, the compiler will repack IF and ONLY IF it needs too.  We may get a performance boost by eliminating manual repacking and using `real, contiguous, pointer` instead. [Reference](https://modelingguru.nasa.gov/servlet/JiveServlet/previewBody/1527-102-1-2631/N1729-4.pdf) page 7.


## Feature requests
* (SOMEWHAT IMPORTANT) Parallel domain decomposition in x2 *and* x3 - this is a big task that is likely to be left aside until I can renew funding.  It's also questionable how useful it is at this point where my typical runs are 32-256 cores (although undoubtedly it may become useful for runs with thousands of cores).  I've found good speedup even dividing the x3 dimension into slabs 2 grid points wide; although that means passing essentially all the grid data around via mpi, the large number of operations per slab means that the effective overhead here is not too much to prevent this from being useful.  However, for simulations that run with GLOW this will massively speed things up...
* (SOMEWHAT IMPORTANT) Periodically updating background neutral atmosphere - should really be done for simulations more than a few hours long but will affect performance
* (EFFICIENCY) Exclusion of null points from field aligned advection, thermal conduction, and source terms - could improve performance
* (INITIAL IMPLEMENTATION COMPLETED) HDF5 file input and output
* (INITIAL IMPLEMENTATION COMPLETED) Option to run the code in a single precision mode - would help with memory limited systems although it's not clear how this would impact numerics (I've never tested my methods in single precision)
* (INITIAL IMPLEMENTATION COMPLETED) Add 3Dtest to ctest
* Add a script that runs a complete sequences of and example from generating ICs, to running the disturbance simulation.
* For local simulations it makes sense to have the simulation able to use a Lagrangian frame of reference to reduce the total number of grid points needed.
* Possibly merge in P. Inchin's EIA changes (with appropriate flags)
* Add an example or options to run with a global grid, in case that is useful for anyone.  This basically requires a special grid generation script and then the simulation needs to be flagged as periodic in x3 (magnetic longitude).  
* Ability to run a dipole grid that encapsulates the magnetic poles


## Plans for adding physics:
These are projects in progress involved GEMINI, you are encouraged to email M. Zettergren for more info if you have interest in using or collaborating on these so that we can efficiently combine efforts and avoid duplicative work.

* Resolved potential solutions - decimate parallel grid down to Farley mapping scale for perp resolution then so the solve on that coarse grid then interpolate back up to original grid.  I've had luck with MUMPS solves in reasonable time up to 300 x 300 x 15 grid points which is probably enough to do something interesting with appropriate periodic and Lagrangian grids (moving at E x B).  
* Diamagnetic drift and perpendicular ambipolar fields - necessary for the smallest scales, e.g. less than 100 m
* Need to add option for true coordinates to be used in the computations of magnetic perturbations (instead of flattened-out spherical)


## Interfaces with other models
* (VERY IMPORTANT) Ability to use GLOW to compute ionization and heating rates, as well as brightnesses of various bands of interest.  This brings up a lot of questions about how GLOW will function on a closed field-line grid; we may need to talk to Stan about this.  


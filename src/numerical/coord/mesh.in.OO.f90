module mesh

!> This is a base type for defining most general characteristics of a curvlinear mesh; extension define specific coordinate systems.
!    However, the idea here is to insulate the numerical parts of the program from that so that they can deal exclusively with the
!    generic metric factors, etc.

use, intrinsic:: iso_fortran_env, only: wp=>real@realbits@

implicit none (type, external)
public

!> curvmesh is an overarching dervied type containing functionality and data that is not specific to individual coordinate systems
!(which are extended types)
type :: curvmesh
  !> SIZE INFORMATION.  Specified and set by base class methods
  integer :: lx1,lx2,lx3,lx2all,lx3all
  !! for program units that may not be able to access module globals
  
  !!> CURVILINEAR VARIABLES AND DIFFS.
  real(wp), dimension(:), allocatable :: x1     ! provided in input file
  real(wp), dimension(:), allocatable :: x1i    ! recomputed by base class once x1 set
  real(wp), dimension(:), pointer :: dx1        ! recomputed by base class from x1
  !! because sub arrays need to be assigned to aliases in calculus module program units
  real(wp), dimension(:), allocatable :: dx1i   ! recomputed by base class from interface locations (themselves recomputed)
  
  real(wp), dimension(:), allocatable :: x2
  real(wp), dimension(:), allocatable :: x2i
  real(wp), dimension(:), pointer :: dx2
  real(wp), dimension(:), allocatable :: dx2i
  
  real(wp), dimension(:), allocatable :: x3
  real(wp), dimension(:), allocatable :: x3i
  real(wp), dimension(:), pointer :: dx3
  real(wp), dimension(:), allocatable :: dx3i
  
  real(wp), dimension(:), allocatable :: x2all
  real(wp), dimension(:), allocatable  :: x2iall
  real(wp), dimension(:), pointer :: dx2all
  real(wp), dimension(:), allocatable  :: dx2iall
  
  real(wp), dimension(:), allocatable :: x3all
  real(wp), dimension(:), allocatable  :: x3iall
  real(wp), dimension(:), pointer :: dx3all
  real(wp), dimension(:), allocatable  :: dx3iall
  
  !> DIFFERENTIAL LENGTH ELEMENTS NEEDED TO COMPUTE COURANT NUMBERS.  Compute by a method in generic class
  real(wp), dimension(:,:,:), allocatable :: dl1i,dl2i,dl3i
  
  !> A FLAG FOR INDICATING WHETHER OR NOT PERIODIC.  set by input files
  logical :: flagper
  
  !> flag for indicated type of grid (0 - closed dipole; 1 - open dipole inverted; 2 - non-inverted).  Computed by method in generic
  !class
  integer :: gridflag
  
  !> METRIC FACTORS.  These are pointers to be assigned/allocated/filled by subclass for specific coordinate system
  real(wp), dimension(:,:,:), pointer :: h1,h2,h3     ! need to be computed by subclass for specific coordinate system
  !! these are the cell-centered metric coefficients
  real(wp), dimension(:,:,:), allocatable :: h1x1i,h2x1i,h3x1i
  !! metric factors at x1 cell interfaces; dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), allocatable :: h1x2i,h2x2i,h3x2i
  !! metric factors at x2 interfaces; dim. 2 has lx2+1
  real(wp), dimension(:,:,:), allocatable :: h1x3i,h2x3i,h3x3i
  !! metric factors at x3 interfaces; dim. 3 has lx3+1
  
  !ROOT ONLY FULL GRID METRIC FACTORS (WORKERS WILL NOT ALLOCATE).  Again these must have subclass methods to provide allocation and filling
  real(wp), dimension(:,:,:), pointer :: h1all,h2all,h3all
  real(wp), dimension(:,:,:), allocatable :: h1x1iall,h2x1iall,h3x1iall
  !! dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), allocatable :: h1x2iall,h2x2iall,h3x2iall
  !! dim. 2 has lx2all+1
  real(wp), dimension(:,:,:), allocatable :: h1x3iall,h2x3iall,h3x3iall
  !! dim. 3 has lx3all+1
  
  !> Problem-depedent geometric terms that can be precomputed, e.g. for advection and elliptic equations
  
  !> UNIT VECTORS - Pointers.  Subclass methods must allocate and assign values to these.  
  real(wp), dimension(:,:,:,:), allocatable :: e1,e2,e3
  !! unit vectors in curvilinear space (cartesian components)
  real(wp), dimension(:,:,:,:), allocatable :: er,etheta,ephi
  !! spherical unit vectors (cartesian components)
  
  !> GEOMAGNETIC GRID DATA
  real(wp), dimension(:,:,:), allocatable :: r,theta,phi
  !! may be used in the interpolation of neutral perturbations
  
  !> FULL-GRID GEOMAGNETIC INFORMATION - USED BY ROOT IN INTERPOLATING ELECTRIC FIELD FILE INPUT.  Pointers.  Subclass methods for each
  !coordinate system must allocate and compute these. 
  real(wp), dimension(:,:,:), allocatable :: rall,thetaall,phiall
  
  !> GEOGRAPHIC DATA; pointers
  real(wp), dimension(:,:,:), allocatable :: glat,glon,alt
  
  !> MAGNETIC FIELD - THIS IS  PART OF THE GRID SINCE THE COORDINATE SYSTEM USED IS BASED ON THE MAGNETIC FIELD.  Pointers
  real(wp), dimension(:,:,:), allocatable :: Bmag
  real(wp), dimension(:,:), allocatable :: I
  
  !> NEED FOR EIA CALCULATIONS.  Pointers
  real(wp), dimension(:,:,:), allocatable :: altall,glonall
  real(wp), dimension(:,:,:), allocatable :: Bmagall
  
  !> DEFINE POINTS TO EXCLUDE FROM NUMERICAL SOLVES?. Pointers
  real(wp), dimension(:,:,:), allocatable :: nullpts
  !! this could be a logical but I'm going to treat it as real*8
  integer :: lnull
  !! length of null point index array
  integer, dimension(:,:), allocatable :: inull

  contains
    procedure :: calc_coord_diffs
    procedure :: set_sizes
    procedure :: set_coords
    procedure :: calc_difflengths
    final :: destructor
end type curvmesh

end module mesh

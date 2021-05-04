module mesh

use, intrinsic:: iso_fortran_env, only: wp=>real@realbits@

implicit none (type, external)
public

type :: curvmesh
!! CURVILINEAR VARIABLES AND DIFFS.
real(wp), dimension(:), allocatable :: x1
real(wp), dimension(:), allocatable :: x1i
real(wp), dimension(:), pointer :: dx1
!! because sub arrays need to be assigned to aliases in calculus module program units
real(wp), dimension(:), allocatable :: dx1i

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

!METRIC FACTORS
real(wp), dimension(:,:,:), pointer :: h1,h2,h3
!! these are the cell-centered metric coefficients
real(wp), dimension(:,:,:), allocatable :: h1x1i,h2x1i,h3x1i
!! metric factors at x1 cell interfaces; dimension 1 has size lx1+1
real(wp), dimension(:,:,:), allocatable :: h1x2i,h2x2i,h3x2i
!! metric factors at x2 interfaces; dim. 2 has lx2+1
real(wp), dimension(:,:,:), allocatable :: h1x3i,h2x3i,h3x3i
!! metric factors at x3 interfaces; dim. 3 has lx3+1

!ROOT ONLY FULL GRID METRIC FACTORS (WORKERS WILL NOT ALLOCATE)
real(wp), dimension(:,:,:), pointer :: h1all,h2all,h3all
real(wp), dimension(:,:,:), allocatable :: h1x1iall,h2x1iall,h3x1iall
!! dimension 1 has size lx1+1
real(wp), dimension(:,:,:), allocatable :: h1x2iall,h2x2iall,h3x2iall
!! dim. 2 has lx2all+1
real(wp), dimension(:,:,:), allocatable :: h1x3iall,h2x3iall,h3x3iall
!! dim. 3 has lx3all+1

!SHALL WE ALSO PRECOMPUTE SOME OF THE PRODUCTS FOR ADVECTION?

!> SIZE INFORMATION
integer :: lx1,lx2,lx3,lx2all,lx3all
!! for program units that may not be able to access module globals

!> UNIT VECTORS - NO ONE NEEDS A FULL GRID COPY OF THESE, I BELIEVE
real(wp), dimension(:,:,:,:), allocatable :: e1,e2,e3
!! unit vectors in curvilinear space (in cartesian components)
real(wp), dimension(:,:,:,:), allocatable :: er,etheta,ephi
!! spherical unit vectors (in cartesian components)

!> GEOMAGNETIC GRID DATA
real(wp), dimension(:,:,:), allocatable :: r,theta,phi
!! may be used in the interpolation of neutral perturbations

!> FULL-GRID GEOMAGNETIC INFORMATION - USED BY ROOT IN INTERPOLATING ELECTRIC FIELD FILE INPUT
real(wp), dimension(:,:,:), allocatable :: rall,thetaall,phiall

!> GEOGRAPHIC DATA
real(wp), dimension(:,:,:), allocatable :: glat,glon,alt

!> MAGNETIC FIELD - THIS IS  PART OF THE GRID SINCE THE COORDINATE SYSTEM USED IS BASED ON THE MAGNETIC FIELD
real(wp), dimension(:,:,:), allocatable :: Bmag
real(wp), dimension(:,:), allocatable :: I

!> NEED FOR EIA CALCULATIONS
real(wp), dimension(:,:,:), allocatable :: altall,glonall
real(wp), dimension(:,:,:), allocatable :: Bmagall

!> DEFINE POINTS TO EXCLUDE FROM NUMERICAL SOLVES?
real(wp), dimension(:,:,:), allocatable :: nullpts
!! this could be a logical but I'm going to treat it as real*8
integer :: lnull
!! length of null point index array
integer, dimension(:,:), allocatable :: inull

!> DIFFERENTIAL LENGTH ELEMENTS NEEDED TO COMPUTE COURANT NUMBERS
real(wp), dimension(:,:,:), allocatable :: dl1i,dl2i,dl3i

!> A FLAG FOR INDICATING WHETHER OR NOT PERIODIC
logical :: flagper

!> flag for indicated type of grid (0 - closed dipole; 1 - open dipole inverted; 2 - non-inverted)
integer :: gridflag

end type curvmesh

end module mesh

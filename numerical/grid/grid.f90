module grid
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

use mpi, only: mpi_integer, mpi_comm_world, mpi_status_ignore

use phys_consts, only: Gconst,Me,Re,wp, red, black

use mpimod, only: myid, lid, &
  tagx1, tagx2, tagx3, tagtheta, tagr, tagphi, tagnull, taglx1, taglx2, taglx3, taglx3all, taginc, &
  tagh1, tagh2, tagh3, tagglat, tagglon, tageunit1, tageunit2, tageunit3, tagetheta, tager, &
  tagalt, tagbmag, tagephi, tagswap, &
  mpi_realprec, ierr, &
  bcast_recv, bcast_send, bcast_recv3d_ghost, bcast_send3d_ghost, bcast_recv3d_x3i, bcast_send3d_x3i

implicit none

private

integer, protected :: lx1,lx2,lx3,lx3all    !this is a useful shorthand for most program units using this module, occassionally a program unit needs to define its own size in which case an only statement is required when using this module.

integer, private :: inunit, &     ! for grid input file
                    outunit       ! for debugging grid (output file for testing)


real(wp), dimension(:,:,:), allocatable, protected :: g1,g2,g3   !gravity, not to be modified by a procedure outside this module
integer, protected :: gridflag    !for cataloguing the type of grid that we are using, open, closed, inverted, etc.
integer :: flagswap    !have the x2 and x3 dimensions been swapped?


type :: curvmesh
  !CURVILINEAR VARIABLES AND DIFFS.
  real(wp), dimension(:), allocatable :: x1
  real(wp), dimension(:), allocatable :: x1i
  real(wp), dimension(:), pointer :: dx1    !because sub arrays need to be assigned to aliases in calculus module program units
  real(wp), dimension(:), allocatable :: dx1i

  real(wp), dimension(:), allocatable :: x2
  real(wp), dimension(:), allocatable :: x2i
  real(wp), dimension(:), pointer :: dx2
  real(wp), dimension(:), allocatable :: dx2i

  real(wp), dimension(:), allocatable :: x3
  real(wp), dimension(:), allocatable :: x3i
  real(wp), dimension(:), pointer :: dx3
  real(wp), dimension(:), allocatable :: dx3i

  real(wp), dimension(:), allocatable :: x3all
  real(wp), dimension(:), allocatable  :: x3iall
  real(wp), dimension(:), pointer :: dx3all
  real(wp), dimension(:), allocatable  :: dx3iall

  !METRIC FACTORS
  real(wp), dimension(:,:,:), pointer :: h1,h2,h3            !these are the cell-centered metric coefficients
  real(wp), dimension(:,:,:), allocatable :: h1x1i,h2x1i,h3x1i   !metric factors at x1 cell interfaces; dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), allocatable :: h1x2i,h2x2i,h3x2i   !metric factors at x2 interfaces; dim. 2 has lx2+1
  real(wp), dimension(:,:,:), allocatable :: h1x3i,h2x3i,h3x3i   !metric factors at x3 interfaces; dim. 3 has lx3+1

  !ROOT ONLY FULL GRID METRIC FACTORS (WORKERS WILL NOT ALLOCATE) - CANDIDATE FOR ELIMINATION???
  real(wp), dimension(:,:,:), pointer :: h1all,h2all,h3all
  real(wp), dimension(:,:,:), allocatable :: h1x1iall,h2x1iall,h3x1iall   !dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), allocatable :: h1x2iall,h2x2iall,h3x2iall   !dim. 2 has lx2+1
  real(wp), dimension(:,:,:), allocatable :: h1x3iall,h2x3iall,h3x3iall   !dim. 3 has lx3all+1

  !SHALL WE ALSO PRECOMPUTE SOME OF THE PRODUCTS FOR ADVECTION?

  !SIZE INFORMATION
  integer :: lx1,lx2,lx3,lx3all    !for program units that may not be able to access module globals

  !UNIT VECTORS - NO ONE NEEDS A FULL GRID COPY OF THESE, I BELIEVE
  real(wp), dimension(:,:,:,:), allocatable :: e1,e2,e3     !unit vectors in curvilinear space (in cartesian components)
  real(wp), dimension(:,:,:,:), allocatable :: er,etheta,ephi    !spherical unit vectors (in cartesian components)

  !GEOMAGNETIC GRID DATA
  real(wp), dimension(:,:,:), allocatable :: r,theta,phi    !may be used in the interpolation of neutral perturbations

  !FULL-GRID GEOMAGNETIC INFORMATION - USED BY ROOT IN INTERPOLATING ELECTRIC FIELD FILE INPUT
  real(wp), dimension(:,:,:), allocatable :: rall,thetaall,phiall

  !GEOGRAPHIC DATA
  real(wp), dimension(:,:,:), allocatable :: glat,glon,alt

  !MAGNETIC FIELD - THIS IS  PART OF THE GRID SINCE THE COORDINATE SYSTEM USED IS BASED ON THE MAGNETIC FIELD
  real(wp), dimension(:,:,:), allocatable :: Bmag
  real(wp), dimension(:,:), allocatable :: I

  !DEFINE POINTS TO EXCLUDE FROM NUMERICAL SOLVES?
  real(wp), dimension(:,:,:), allocatable :: nullpts   !this could be a logical but I'm going to treat it as real*8
  integer :: lnull                                    !length of null point index array
  integer, dimension(:,:), allocatable :: inull

  !DIFFERENTIAL LENGTH ELEMENTS NEEDED TO COMPUTE COURANT NUMBERS
  real(wp), dimension(:,:,:), allocatable :: dl1i,dl2i,dl3i

  !A FLAG FOR INDICATING WHETHER OR NOT PERIODIC
  integer :: flagper
end type curvmesh


public :: curvmesh,  lx1,lx2,lx3, lx3all, gridflag, flagswap, clear_unitvecs, g1,g2,g3, &
  read_grid, clear_grid

contains


subroutine read_grid(indatsize,indatgrid,flagperiodic,x)

!------------------------------------------------------------
!--------READS GRID INFORMATION FROM A BBINARY FILE AND ALLOCATES
!--------GRID STRUCTURE.  NOTE THAT THE INPUT FILE DOES NOT,
!--------BY DEFAULT INCLUDE THE GHOST CELLS.  THIS FUNCTION
!--------IS A WRAPPER WHICH PASSES WORK OFF TO APPROPRIATE
!--------SUBROUTINE DEPENDING ON WHETHER IT IS CALLED BY ROOT
!--------OR NOT.   THIS SUBRTOUINE ALSO CLASSIFIES THE GRID.
!------------------------------------------------------------

character(*), intent(in) :: indatsize,indatgrid
integer, intent(in) :: flagperiodic
type(curvmesh), intent(inout) :: x


!READ IN THE GRID DATA
if(myid==0) then
  call read_grid_root(indatsize,indatgrid,x)
else
  call read_grid_workers(x)
end if


!FLAG THE GRID AS PERIODIC, IF REQUESTED; IF SO PERIODICITY WILL BE ASSUMED
!IN THE X3-DIRECTION, NOTE BOTH ROOT AND WORKERS MUST DO THIS!!!
if (flagperiodic==1) then
  x%flagper=1
else
  x%flagper=0
end if


!DETERMINE THE TYPE OF GRID WE HAVE AND SET AN APPROPRIATE FLAG
if (abs(x%alt(1,1,1)-x%alt(lx1,1,1))<10d3) then    !closed dipole grid
  gridflag=0
else if (x%alt(1,1,1)>x%alt(2,1,1)) then    !open dipole grid with inverted structure wrt altitude
  gridflag=1
else    !something different (viz. non-inverted - lowest altitudes at the logical bottom of the grid)
  gridflag=2
end if

end subroutine read_grid


subroutine read_grid_root(indatsize,indatgrid,x)

!------------------------------------------------------------
!--------SOME CODE DUPLICATION WITH WORKER VERSION - CAN WE
!--------CREATE A COMMON SUBROUTINE THAT ALLOCATES SOME VARS?
!------------------------------------------------------------

character(*), intent(in) :: indatsize,indatgrid
type(curvmesh), intent(inout) :: x    !does this need to be inout?  I think if anything in unallocated, it does...

integer lx1g,lx2g,lx3allg,iid,ix1,ix2,ix3,icount,icomp!,itell

!NOTE THAT HAVING THESE ARE LOCAL (TEMPORARY) VARS. PREVENTS ROOT FROM WRITING ENTIRE GRID TO FILE AT SOME LATER POINT...
real(wp), dimension(:,:,:), allocatable :: g1all,g2all,g3all   !to temporarily store input data to be distributed
real(wp), dimension(:,:,:), allocatable :: altall,glatall,glonall
real(wp), dimension(:,:,:), allocatable :: rall,thetaall,phiall
real(wp), dimension(:,:,:), allocatable :: Bmagall
real(wp), dimension(:,:), allocatable :: Incall
real(wp), dimension(:,:,:), allocatable :: nullptsall
real(wp), dimension(:,:,:,:), allocatable :: e1all,e2all,e3all,erall,ethetaall,ephiall    !might be best to have a tmp vector array...
real(wp), dimension(:,:,:), allocatable :: mpisendbuf
real(wp), dimension(:,:,:), allocatable :: mpirecvbuf
real(wp), dimension(:,:,:), allocatable :: tmpdx
real(wp), dimension(:,:,:), allocatable :: htmp
real(wp), dimension(:,:), allocatable :: htmp2D
real(wp), dimension(:,:,:,:), allocatable :: htmp4D

logical exists

inquire(file=indatsize, exist=exists)
if (.not.exists) then
   write(stderr,*) 'must generate grid with script before running simulation--grid not present: ',indatsize
   error stop 
endif

!DETERMINE THE SIZE OF THE GRID TO BE LOADED
open(newunit=inunit,file=indatsize,status='old',form='unformatted', &
     access='stream', action='read')
read(inunit) lx1g,lx2g,lx3allg    !note that these are sizes *including ghost cells*  !MZ - lx2allg+declaration above
close(inunit)


!MZ must set up lx2all and lx2 here
!DETERMINE NUMBER OF SLABS AND CORRESPONDING SIZE FOR EACH WORKER 
!NOTE THAT WE WILL ASSUME THAT THE GRID SIZE IS DIVISIBLE BY NUMBER OF WORKERS AS THIS HAS BEEN CHECKED A FEW LINES BELOW
x%lx1=lx1g; x%lx2=lx2g; x%lx3all=lx3allg;
lx1=x%lx1; lx2=x%lx2; lx3all=x%lx3all;    !define some shortand variables for ease of reference
lx3=lx3all/lid                            !note the integer division here   !MZ-change to lid3
x%lx3=lx3
print *, 'Grid has size:  ',lx1,lx2,lx3all    !MZ - lx2all
print *, '    slab size:  ',lx1,lx2,lx3


!ADJUST THE SIZES OF THE VARIABLES IF LX3ALL==1, SO THAT THE ALLOCATIONS ARE THE CORRECT SIZE
if (lx3all==1) then
  print *, 'Detected a 2D run request...  swapping x2 an x3 sizes to maintain parallelization.'
  lx3all=lx2; x%lx3all=lx2;
  lx2=1; x%lx2=1;
  lx3=lx3all/lid
  x%lx3=lx3
  flagswap=1
else
  if(lx2==1) then
    print *, 'Detected a 2D run request...  ***NOT*** swapping x2 an x3 sizes to maintain parallelization.'
  else
    print *, 'Detected a 3D run requuest...'
  endif
  flagswap=0
end if


!JUST BAIL ON THE SIMULATION IF THE X3 SIZE ISN'T VISIBLE BY NUMBER OF PROCESSES
if (lx3*lid/=lx3all) then
  write(stderr,*) 'Number of grid points', lx3all, &
    'must be divisible by number of processes',lid
  error stop 
end if
if (lx3<2) then
  write(stderr,*) red // '**************************************************************************'
  write(stderr,*) 'WARNING/ERROR: simulation with slab size < 2 may give incorrect results with some MPI versions. '
  write(stderr,*) 'Check results on system with MPI -np >= 2.   here, lx3=',lx3
  write(stderr,*) '**************************************************************************' // black
  error stop
end if

!COMMUNICATE THE GRID SIZE TO THE WORKERS SO THAT THEY CAN ALLOCATE SPACE
do iid=1,lid-1
  call mpi_send(lx1,1,MPI_INTEGER,iid,taglx1,MPI_COMM_WORLD,ierr)
  call mpi_send(lx2,1,MPI_INTEGER,iid,taglx2,MPI_COMM_WORLD,ierr)          !need to also pass the lx2all size to all workers to they know
  call mpi_send(lx3,1,MPI_INTEGER,iid,taglx3,MPI_COMM_WORLD,ierr)
  call mpi_send(lx3all,1,MPI_INTEGER,iid,taglx3all,MPI_COMM_WORLD,ierr)    !not clear whether workers actually need the x3all variable, hopefully not, but they definitely need to know the full grid size for the derivative functions
end do


!TELL WORKERS IF WE'VE SWAPPED DIMENSIONS
do iid=1,lid-1
  call mpi_send(flagswap,1,MPI_INTEGER,iid,tagswap,MPI_COMM_WORLD,ierr)
end do


!ALLOCATE SPACE FOR ROOTS SLAB OF DATA
allocate(x%x1(-1:lx1+2))
allocate(x%dx1i(lx1),x%x1i(lx1+1),x%dx1(0:lx1+2))
allocate(x%x2(-1:lx2+2))    !MZ - x2 needs to be broken out like x3all immediately below
allocate(x%dx2i(lx2),x%x2i(lx2+1),x%dx2(0:lx2+2))  !MZ - moved below with dx3, etc.


!FULL-GRID X3-VARIABLE
allocate(x%x3all(-1:lx3all+2))
allocate(x%dx3all(0:lx3all+2))
allocate(x%x3iall(lx3all+1),x%dx3iall(lx3all))

allocate(x%h1all(-1:lx1+2,-1:lx2+2,-1:lx3all+2),x%h2all(-1:lx1+2,-1:lx2+2,-1:lx3all+2), &
         x%h3all(-1:lx1+2,-1:lx2+2,-1:lx3all+2))    !do we need the ghost cell values?  Yes the divergence in compression term is computed using one ghost cell in each direction!!!!
allocate(x%h1x1iall(1:lx1+1,1:lx2,1:lx3all),x%h2x1iall(1:lx1+1,1:lx2,1:lx3all),x%h3x1iall(1:lx1+1,1:lx2,1:lx3all))
allocate(x%h1x2iall(1:lx1,1:lx2+1,1:lx3all),x%h2x2iall(1:lx1,1:lx2+1,1:lx3all),x%h3x2iall(1:lx1,1:lx2+1,1:lx3all))
allocate(x%h1x3iall(1:lx1,1:lx2,1:lx3all+1),x%h2x3iall(1:lx1,1:lx2,1:lx3all+1),x%h3x3iall(1:lx1,1:lx2,1:lx3all+1))

allocate(x%rall(1:lx1,1:lx2,1:lx3all),x%thetaall(1:lx1,1:lx2,1:lx3all),x%phiall(1:lx1,1:lx2,1:lx3all))


!DETERMINE AND ALLOCATE SPACE NEEDED FOR ROOT SUBGRIDS (WORKERS USE SIMILAR ALLOCATE STATEMENTS)
allocate(x%x3(-1:lx3+2))
allocate(x%dx3i(lx3),x%x3i(lx3+1),x%dx3(0:lx3+2))

allocate(x%h1(-1:lx1+2,-1:lx2+2,-1:lx3+2),x%h2(-1:lx1+2,-1:lx2+2,-1:lx3+2),x%h3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
!    allocate(x%h1(1:lx1,1:lx2,1:lx3),x%h2(1:lx1,1:lx2,1:lx3),x%h3(1:lx1,1:lx2,1:lx3))
allocate(x%h1x1i(1:lx1+1,1:lx2,1:lx3),x%h2x1i(1:lx1+1,1:lx2,1:lx3),x%h3x1i(1:lx1+1,1:lx2,1:lx3))
allocate(x%h1x2i(1:lx1,1:lx2+1,1:lx3),x%h2x2i(1:lx1,1:lx2+1,1:lx3),x%h3x2i(1:lx1,1:lx2+1,1:lx3))
allocate(x%h1x3i(1:lx1,1:lx2,1:lx3+1),x%h2x3i(1:lx1,1:lx2,1:lx3+1),x%h3x3i(1:lx1,1:lx2,1:lx3+1))

allocate(x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3),x%alt(1:lx1,1:lx2,1:lx3))
allocate(x%r(1:lx1,1:lx2,1:lx3),x%theta(1:lx1,1:lx2,1:lx3),x%phi(1:lx1,1:lx2,1:lx3))

allocate(x%Bmag(1:lx1,1:lx2,1:lx3))
allocate(x%I(1:lx2,1:lx3))
allocate(x%nullpts(1:lx1,1:lx2,1:lx3))

allocate(x%e1(1:lx1,1:lx2,1:lx3,1:3),x%e2(1:lx1,1:lx2,1:lx3,1:3),x%e3(1:lx1,1:lx2,1:lx3,1:3))
allocate(x%er(1:lx1,1:lx2,1:lx3,1:3),x%etheta(1:lx1,1:lx2,1:lx3,1:3),x%ephi(1:lx1,1:lx2,1:lx3,1:3))


!NOW WE NEED TO READ IN THE FULL GRID DATA AND PUT IT INTO THE STRUCTURE.
!IF WE HAVE DONE ANY DIMENSION SWAPPING HERE WE NEED TO TAKE THAT INTO ACCOUNT IN THE VARIABLES THAT ARE BEING READ IN
print *, 'Starting grid input from file: ',indatgrid
open(newunit=inunit,file=indatgrid,status='old',form='unformatted',access='stream', action='read')
if (flagswap/=1) then                     !normal (i.e. full 3D) grid ordering, or a 2D grid with 1 element naturally in the second dimension
  read(inunit) x%x1,x%x1i,x%dx1,x%dx1i
  read(inunit) x%x2,x%x2i,x%dx2,x%dx2i     !MZ - these becomes 'all' variables
  read(inunit) x%x3all,x%x3iall,x%dx3all,x%dx3iall
  read(inunit) x%h1all,x%h2all,x%h3all
  read(inunit) x%h1x1iall,x%h2x1iall,x%h3x1iall
  read(inunit) x%h1x2iall,x%h2x2iall,x%h3x2iall
  read(inunit) x%h1x3iall,x%h2x3iall,x%h3x3iall

  allocate(g1all(lx1,lx2,lx3all),g2all(lx1,lx2,lx3all),g3all(lx1,lx2,lx3all))    !MZ - lx2all
  read(inunit) g1all,g2all,g3all

  allocate(altall(lx1,lx2,lx3all),glatall(lx1,lx2,lx3all),glonall(lx1,lx2,lx3all))    !MZ - lx2all
  read(inunit) altall,glatall,glonall

  allocate(Bmagall(lx1,lx2,lx3all))     !MZ - lx2all
  read(inunit) Bmagall

  allocate(Incall(lx2,lx3all))          !MZ - lx2all
  read(inunit) Incall

  allocate(nullptsall(lx1,lx2,lx3all))   !MZ - lx2all
  read(inunit) nullptsall


  !MZ - remainder require lx2all
  allocate(e1all(lx1,lx2,lx3all,3))
  read(inunit) e1all
  allocate(e2all(lx1,lx2,lx3all,3))
  read(inunit) e2all
  allocate(e3all(lx1,lx2,lx3all,3))
  read(inunit) e3all

  allocate(erall(lx1,lx2,lx3all,3))
  read(inunit) erall
  allocate(ethetaall(lx1,lx2,lx3all,3))
  read(inunit) ethetaall
  allocate(ephiall(lx1,lx2,lx3all,3))
  read(inunit) ephiall

  allocate(rall(lx1,lx2,lx3all))
  read(inunit) rall
  allocate(thetaall(lx1,lx2,lx3all))
  read(inunit) thetaall
  allocate(phiall(lx1,lx2,lx3all))
  read(inunit) phiall

  x%rall=rall; x%thetaall=thetaall; x%phiall=phiall;
else     !this is apparently a 2D grid, so the x2 and x3 dimensions have been/need to be swapped   !MZ - may need to change lx2-->lx2all???
  print *, 'Detected a 2D grid, so will permute the dimensions of the input'

  !Note there that the fortran arrays are the correct size, but the input data are not!!!  This means tmp variable and permutes...
  read(inunit) x%x1,x%x1i,x%dx1,x%dx1i                !x1 untouched
  read(inunit) x%x3all,x%x3iall,x%dx3all,x%dx3iall    !for a 3D grid this is x2, but now considered x3(all)
  read(inunit) x%x2,x%x2i,x%dx2,x%dx2i                !formerly x3, now x2

  allocate(htmp(-1:lx1+2,-1:lx3all+2,-1:lx2+2))    !this stores the input metric factors which are swapped x2/x3 vs. what this simulation will use
  read(inunit) htmp
  x%h1all=reshape(htmp,[lx1+4,lx2+4,lx3all+4],order=[1,3,2])
  read(inunit) htmp        !this would be h3, but with the input structure shape
  x%h3all=reshape(htmp,[lx1+4,lx2+4,lx3all+4],order=[1,3,2])   !permute the dimensions of the array 3 --> 2, 2 --> 3
  read(inunit) htmp        !this would be h3, but with the input structure shape
  x%h2all=reshape(htmp,[lx1+4,lx2+4,lx3all+4],order=[1,3,2])
  deallocate(htmp)

  allocate(htmp(1:lx1+1,1:lx3all,1:lx2))    !input 2 vs. 3 dimensions swapped from this program
  read(inunit) htmp
  x%h1x1iall=reshape(htmp,[lx1+1,lx2,lx3all],order=[1,3,2])
  read(inunit) htmp
  x%h3x1iall=reshape(htmp,[lx1+1,lx2,lx3all],order=[1,3,2])
  read(inunit) htmp
  x%h2x1iall=reshape(htmp,[lx1+1,lx2,lx3all],order=[1,3,2])
  deallocate(htmp)

  allocate(htmp(1:lx1,1:lx3all+1,1:lx2))
  read(inunit) htmp
  x%h1x3iall=reshape(htmp,[lx1,lx2,lx3all+1],order=[1,3,2])    !Note also that the x2 interface from teh input file is x3i in this simulation
  read(inunit) htmp
  x%h3x3iall=reshape(htmp,[lx1,lx2,lx3all+1],order=[1,3,2])
  read(inunit) htmp
  x%h2x3iall=reshape(htmp,[lx1,lx2,lx3all+1],order=[1,3,2])
  deallocate(htmp)

  allocate(htmp(1:lx1,1:lx3all,1:lx2+1))
  read(inunit) htmp
  x%h1x2iall=reshape(htmp,[lx1,lx2+1,lx3all],order=[1,3,2])
  read(inunit) htmp
  x%h3x2iall=reshape(htmp,[lx1,lx2+1,lx3all],order=[1,3,2])
  read(inunit) htmp
  x%h2x2iall=reshape(htmp,[lx1,lx2+1,lx3all],order=[1,3,2])
  deallocate(htmp)

  allocate(g1all(lx1,lx2,lx3all),g2all(lx1,lx2,lx3all),g3all(lx1,lx2,lx3all))
  allocate(htmp(lx1,lx3all,lx2))
  read(inunit) htmp
  g1all=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  read(inunit) htmp
  g3all=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  read(inunit) htmp
  g2all=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  deallocate(htmp)

  allocate(altall(lx1,lx2,lx3all),glatall(lx1,lx2,lx3all),glonall(lx1,lx2,lx3all))
  allocate(htmp(lx1,lx3all,lx2))
  read(inunit) htmp
  altall=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])      
  read(inunit) htmp
  glatall=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  read(inunit) htmp
  glonall=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  deallocate(htmp)

  allocate(Bmagall(lx1,lx2,lx3all))
  allocate(htmp(lx1,lx3all,lx2))
  read(inunit) htmp
  Bmagall=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  deallocate(htmp)

  allocate(Incall(lx2,lx3all))
  allocate(htmp2D(lx3all,lx2))
  read(inunit) htmp2D
  Incall=reshape(htmp2D,[lx2,lx3all],order=[2,1])
  deallocate(htmp2D)


  !inquire(inunit, pos=itell)
  !print *,'file pos before read',itell
  !read(inunit) nullptsall
  !print *,'lx1',lx1,'lx2',lx2,'lx3all',lx3all
  !print *,'shape(nullptsall)',shape(nullptsall)
  !inquire(inunit, pos=itell)
  !print *,'file pos after read',itell
! FIXME BROKEN!
  !allocate(nullptsall(lx1,lx2,lx3all))
  !print *,shape(nullptsall)
  !stop
  
  ! FIXME would be like this, but this doesn't work.
  !allocate(nullptsall(lx1,lx2,lx3all))
  !read(inunit) nullptsall

  allocate(htmp(lx1,lx3all,lx2))
  read(inunit) htmp
  nullptsall=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  deallocate(htmp)

  allocate(e1all(lx1,lx2,lx3all,3))
  allocate(htmp4D(lx1,lx3all,lx2,3))
  read(inunit) htmp4D
  e1all=reshape(htmp4D,[lx1,lx2,lx3all,3],order=[1,3,2,4])
  deallocate(htmp4D)

  allocate(e3all(lx1,lx2,lx3all,3))          !swap the x2/x3 unit vectors
  allocate(htmp4D(lx1,lx3all,lx2,3))
  read(inunit) htmp4D
  e3all=reshape(htmp4D,[lx1,lx2,lx3all,3],order=[1,3,2,4])
  deallocate(htmp4D)

  allocate(e2all(lx1,lx2,lx3all,3))
  allocate(htmp4D(lx1,lx3all,lx2,3))
  read(inunit) htmp4D
  e2all=reshape(htmp4D,[lx1,lx2,lx3all,3],order=[1,3,2,4])
  deallocate(htmp4D)

  allocate(erall(lx1,lx2,lx3all,3))
  allocate(htmp4D(lx1,lx3all,lx2,3))
  read(inunit) htmp4D
  erall=reshape(htmp4D,[lx1,lx2,lx3all,3],order=[1,3,2,4])
  deallocate(htmp4D)

  allocate(ethetaall(lx1,lx2,lx3all,3))
  allocate(htmp4D(lx1,lx3all,lx2,3))
  read(inunit) htmp4D
  ethetaall=reshape(htmp4D,[lx1,lx2,lx3all,3],order=[1,3,2,4])
  deallocate(htmp4D)

  allocate(ephiall(lx1,lx2,lx3all,3))
  allocate(htmp4D(lx1,lx3all,lx2,3))
  read(inunit) htmp4D
  ephiall=reshape(htmp4D,[lx1,lx2,lx3all,3],order=[1,3,2,4])
  deallocate(htmp4D)

  allocate(rall(lx1,lx2,lx3all),thetaall(lx1,lx2,lx3all),phiall(lx1,lx2,lx3all))
  allocate(htmp(lx1,lx3all,lx2))
  read(inunit) htmp
  rall=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  read(inunit) htmp
  thetaall=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  read(inunit) htmp
  phiall=reshape(htmp,[lx1,lx2,lx3all],order=[1,3,2])
  deallocate(htmp)

  x%rall=rall; x%thetaall=thetaall; x%phiall=phiall;
end if
close(inunit)
print *, 'Done reading in grid data...'


!ALLOCATE SPACE FOR ROOTS SUBGRID GRAVITATIONAL FIELD
allocate(g1(1:lx1,1:lx2,1:lx3),g2(1:lx1,1:lx2,1:lx3),g3(1:lx1,1:lx2,1:lx3))


!SEND FULL X1 AND X2 GRIDS TO EACH WORKER (ONLY X3-DIM. IS INVOLVED IN THE MPI
do iid=1,lid-1
  call mpi_send(x%x1,lx1+4,mpi_realprec,iid,tagx1,MPI_COMM_WORLD,ierr)
  call mpi_send(x%x2,lx2+4,mpi_realprec,iid,tagx2,MPI_COMM_WORLD,ierr)          !MZ - needs to be x2all
  call mpi_send(x%x3all,lx3all+4,mpi_realprec,iid,tagx3,MPI_COMM_WORLD,ierr)    !workers may need a copy of this, e.g. for boudnary conditiosn
  call mpi_send(x%dx1,lx1+3,mpi_realprec,iid,tagx1,MPI_COMM_WORLD,ierr)
  call mpi_send(x%dx2,lx2+3,mpi_realprec,iid,tagx2,MPI_COMM_WORLD,ierr)
  call mpi_send(x%x1i,lx1+1,mpi_realprec,iid,tagx1,MPI_COMM_WORLD,ierr)
  call mpi_send(x%x2i,lx2+1,mpi_realprec,iid,tagx2,MPI_COMM_WORLD,ierr)         !MZ - needs to be all
  call mpi_send(x%dx1i,lx1,mpi_realprec,iid,tagx1,MPI_COMM_WORLD,ierr)
  call mpi_send(x%dx2i,lx2,mpi_realprec,iid,tagx2,MPI_COMM_WORLD,ierr)          !MZ - 'all'
end do
print *, 'Done sending common variables to workers...'


!NOW SEND THE INFO THAT DEPENDS ON X3 SLAB SIZE
call bcast_send(x%x3all,tagx3,x%x3)

x%dx3=x%x3(0:lx3+2)-x%x3(-1:lx3+1)     !computing these avoids extra message passing (could be done for other coordinates, as well)
x%x3i(1:lx3+1)=0.5*(x%x3(0:lx3)+x%x3(1:lx3+1))
x%dx3i=x%x3i(2:lx3+1)-x%x3i(1:lx3)

!MZ - need to have x2 direction dealt with here.  

allocate(mpisendbuf(-1:lx1+2,-1:lx2+2,-1:lx3all+2),mpirecvbuf(-1:lx1+2,-1:lx2+2,-1:lx3+2))

mpisendbuf=x%h1all    !since metric factors are pointers they are not gauranteed to be contiguous in memory so pack them into a buffer that is...
call bcast_send3D_ghost(mpisendbuf,tagh1,mpirecvbuf)    !special broadcast subroutine to handle 3D arrays with ghost cells
x%h1=mpirecvbuf       !store roots slab of metric factors in its grid structure
mpisendbuf=x%h2all
call bcast_send3D_ghost(mpisendbuf,tagh2,mpirecvbuf)
x%h2=mpirecvbuf
mpisendbuf=x%h3all
call bcast_send3D_ghost(mpisendbuf,tagh3,mpirecvbuf)
x%h3=mpirecvbuf

deallocate(mpisendbuf,mpirecvbuf)    !we need different sized buffers below

call bcast_send(x%h1x1iall,tagh1,x%h1x1i)    !do the weird sizes here (ie. lx1+1) give issues with MPI?  probably...  No because bcast reads the size off of the variable...
call bcast_send(x%h2x1iall,tagh2,x%h2x1i)
call bcast_send(x%h3x1iall,tagh3,x%h3x1i)

call bcast_send(x%h1x2iall,tagh1,x%h1x2i)    !MZ - these now become interface sends, x2i, and need correspoding mpi routines written
call bcast_send(x%h2x2iall,tagh2,x%h2x2i)
call bcast_send(x%h3x2iall,tagh3,x%h3x2i)

call bcast_send3D_x3i(x%h1x3iall,tagh1,x%h1x3i)
call bcast_send3D_x3i(x%h2x3iall,tagh2,x%h2x3i)
call bcast_send3D_x3i(x%h3x3iall,tagh3,x%h3x3i)

call bcast_send(g1all,tagh1,g1)
call bcast_send(g2all,tagh2,g2)
call bcast_send(g3all,tagh3,g3)

call bcast_send(glatall,tagglat,x%glat)
call bcast_send(glonall,tagglon,x%glon)
call bcast_send(altall,tagalt,x%alt)

call bcast_send(Bmagall,tagBmag,x%Bmag)
call bcast_send(Incall,taginc,x%I)
call bcast_send(nullptsall,tagnull,x%nullpts)

allocate(mpisendbuf(1:lx1,1:lx2,1:lx3all),mpirecvbuf(1:lx1,1:lx2,1:lx3))    !why is buffering used/needed here???
do icomp=1,3
  mpisendbuf=e1all(:,:,:,icomp) 
  call bcast_send(mpisendbuf,tageunit1,mpirecvbuf)
  x%e1(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=e2all(:,:,:,icomp)
  call bcast_send(mpisendbuf,tageunit2,mpirecvbuf)
  x%e2(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=e3all(:,:,:,icomp)
  call bcast_send(mpisendbuf,tageunit3,mpirecvbuf)
  x%e3(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=erall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tager,mpirecvbuf)
  x%er(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=ethetaall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tagetheta,mpirecvbuf)
  x%etheta(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=ephiall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tagephi,mpirecvbuf)
  x%ephi(:,:,:,icomp)=mpirecvbuf
end do
deallocate(mpisendbuf,mpirecvbuf)

call bcast_send(rall,tagr,x%r)
call bcast_send(thetaall,tagtheta,x%theta)
call bcast_send(phiall,tagphi,x%phi)
print *,  'Done sending slabbed variables to workers...'


!COUNT THE NUMBER OF NULL GRID POINTS AND GENERATE A LIST OF NULL INDICES FOR LATER USE
x%lnull=0;
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3)>0.5d0) then
        x%lnull=x%lnull+1
      end if
    end do
  end do
end do
allocate(x%inull(1:x%lnull,1:3))

icount=1
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3)>0.5d0) then
        x%inull(icount,:)=[ix1,ix2,ix3]
        icount=icount+1
      end if
    end do
  end do
end do
print *, 'Done computing null grid points...  Process:  ',myid,' has:  ',x%lnull


!COMPUTE DIFFERENTIAL DISTANCES ALONG EACH DIRECTION (TO BE USED IN TIME STEP DETERMINATION...
allocate(x%dl1i(1:lx1,1:lx2,1:lx3),x%dl2i(1:lx1,1:lx2,1:lx3),x%dl3i(1:lx1,1:lx2,1:lx3),tmpdx(1:lx1,1:lx2,1:lx3)) 

tmpdx=spread(spread(x%dx1i,2,lx2),3,lx3)
x%dl1i=tmpdx*x%h1(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx2i,1,lx1),3,lx3)
x%dl2i=tmpdx*x%h2(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx3i,1,lx1),2,lx2)
x%dl3i=tmpdx*x%h3(1:lx1,1:lx2,1:lx3)


!    !THIS CODE BLOCK HAS ROOT "PARROT" THE GRID TO A FILE FOR DEBUGGING PURPOSES.  
!    !THIS BLOCK MUST BE HERE DUE TO THE FACT MANY OF HTE ALL-GRID VARIABLES ARE ONLY
!    !IN SCOPE INSIDE THIS FUNCITON BEFORE THE DEALLOCATE STATEMENT BELOW
!    open(newunit=outunit,file='testsize.dat',status='replace',form='unformatted', &
!         access='stream', action='write')
!    write(outunit) lx1,lx2,lx3all    !note that these sizes exxclude ghost cells
!    close(outunit)
!
!    open(outunit,file='testgrid.dat',status='replace',form='unformatted', &
!         access='stream')
!    write(outunit) x%x1,x%x1i,x%dx1,x%dx1i    !I guess the number of elements is obtained from the size of the arrays??? 
!    write(outunit) x%x2,x%x2i,x%dx2,x%dx2i
!    write(outunit) x%x3all,x%x3iall,x%dx3all,x%dx3iall
!    write(outunit) x%h1all,x%h2all,x%h3all
!    write(outunit) x%h1x1iall,x%h2x1iall,x%h3x1iall
!    write(outunit) x%h1x2iall,x%h2x2iall,x%h3x2iall
!    write(outunit) x%h1x3iall,x%h2x3iall,x%h3x3iall
!    write(outunit) g1all,g2all,g3all
!    write(outunit) altall,glatall,glonall
!    write(outunit) Bmagall
!    write(outunit) Incall
!    write(outunit) nullptsall
!    close(outunit)
!    print *, 'Done creating copy of grid...'
!    

!
!    !FLAG THE GRID AS PERIODIC, IF REQUESTED; IF SO PERIODICITY WILL BE ASSUMED IN THE X3-DIRECTION
!    if (flagperiodic==1) then
!      x%flagper=1
!    else
!      x%flagper=0
!    end if
!

!DEALLOCATE ANY FULL GRID VARIABLES THAT ARE NO LONGER NEEDED
deallocate(g1all,g2all,g3all,altall,glatall,glonall,Bmagall,Incall,nullptsall,tmpdx,rall,thetaall,phiall)

end subroutine read_grid_root


subroutine read_grid_workers(x)

type(curvmesh), intent(inout) :: x

integer :: ix1,ix2,ix3,icount,icomp
real(wp), dimension(:,:,:), allocatable :: mpirecvbuf
real(wp), dimension(:,:,:), allocatable :: tmpdx


!GET ROOTS MESSAGE WITH THE SIZE OF THE GRID WE ARE TO RECEIVE
call mpi_recv(lx1,1,MPI_INTEGER,0,taglx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(lx2,1,MPI_INTEGER,0,taglx2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(lx3,1,MPI_INTEGER,0,taglx3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(lx3all,1,MPI_INTEGER,0,taglx3all,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
x%lx1=lx1; x%lx2=lx2; x%lx3all=lx3all;
x%lx3=lx3


!ROOT NEEDS TO TELL US WHETHER WE'VE SWAPPED DIMENSIONS SINCE THIS AFFECTS HOW CURRENTS ARE COMPUTED
call mpi_recv(flagswap,1,MPI_INTEGER,0,tagswap,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)


!ALLOCATE SPACE FOR MY SLAB OF DATA
allocate(x%x1(-1:lx1+2))
allocate(x%dx1i(lx1),x%x1i(lx1+1),x%dx1(0:lx1+2))
allocate(x%x2(-1:lx2+2))
allocate(x%dx2i(lx2),x%x2i(lx2+1),x%dx2(0:lx2+2))    !MZ - need an x2all variable like x3 is done below

!DETERMINE AND ALLOCATE SPACE NEEDED FOR WORKERS SUBGRIDS
allocate(x%x3(-1:lx3+2))
allocate(x%dx3i(lx3),x%x3i(lx3+1),x%dx3(0:lx3+2))
allocate(x%x3all(-1:lx3all+2))

allocate(x%h1(-1:lx1+2,-1:lx2+2,-1:lx3+2),x%h2(-1:lx1+2,-1:lx2+2,-1:lx3+2),x%h3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
allocate(x%h1x1i(1:lx1+1,1:lx2,1:lx3),x%h2x1i(1:lx1+1,1:lx2,1:lx3),x%h3x1i(1:lx1+1,1:lx2,1:lx3))
allocate(x%h1x2i(1:lx1,1:lx2+1,1:lx3),x%h2x2i(1:lx1,1:lx2+1,1:lx3),x%h3x2i(1:lx1,1:lx2+1,1:lx3))
allocate(x%h1x3i(1:lx1,1:lx2,1:lx3+1),x%h2x3i(1:lx1,1:lx2,1:lx3+1),x%h3x3i(1:lx1,1:lx2,1:lx3+1))

allocate(x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3),x%alt(1:lx1,1:lx2,1:lx3))
allocate(x%r(1:lx1,1:lx2,1:lx3),x%theta(1:lx1,1:lx2,1:lx3),x%phi(1:lx1,1:lx2,1:lx3))

allocate(x%Bmag(1:lx1,1:lx2,1:lx3))
allocate(x%I(1:lx2,1:lx3))
allocate(x%nullpts(1:lx1,1:lx2,1:lx3))

allocate(x%e1(1:lx1,1:lx2,1:lx3,1:3),x%e2(1:lx1,1:lx2,1:lx3,1:3),x%e3(1:lx1,1:lx2,1:lx3,1:3))
allocate(x%er(1:lx1,1:lx2,1:lx3,1:3),x%etheta(1:lx1,1:lx2,1:lx3,1:3),x%ephi(1:lx1,1:lx2,1:lx3,1:3))


!ALLOCATE SPACE FOR WORKER'S GRAVITATIONAL FIELD
allocate(g1(1:lx1,1:lx2,1:lx3),g2(1:lx1,1:lx2,1:lx3),g3(1:lx1,1:lx2,1:lx3))


!RECEIVE GRID DATA FROM ROOT
call mpi_recv(x%x1,lx1+4,mpi_realprec,0,tagx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(x%x2,lx2+4,mpi_realprec,0,tagx2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)    !MZ - becomes x2all
call mpi_recv(x%x3all,lx3all+4,mpi_realprec,0,tagx3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(x%dx1,lx1+3,mpi_realprec,0,tagx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(x%dx2,lx2+3,mpi_realprec,0,tagx2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(x%x1i,lx1+1,mpi_realprec,0,tagx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(x%x2i,lx2+1,mpi_realprec,0,tagx2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(x%dx1i,lx1,mpi_realprec,0,tagx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
call mpi_recv(x%dx2i,lx2,mpi_realprec,0,tagx2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

call bcast_recv(x%x3,tagx3)

x%dx3=x%x3(0:lx3+2)-x%x3(-1:lx3+1)     !computing these avoids extra message passing (could be done for other coordinates, as well)
x%x3i(1:lx3+1)=0.5*(x%x3(0:lx3)+x%x3(1:lx3+1))
x%dx3i=x%x3i(2:lx3+1)-x%x3i(1:lx3)

!MZ - need to receive our piece of x2 and computd differences as done for x3 above

allocate(mpirecvbuf(-1:lx1+2,-1:lx2+2,-1:lx3+2))

call bcast_recv3D_ghost(mpirecvbuf,tagh1)
x%h1=mpirecvbuf
call bcast_recv3D_ghost(mpirecvbuf,tagh2)
x%h2=mpirecvbuf
call bcast_recv3D_ghost(mpirecvbuf,tagh3)
x%h3=mpirecvbuf

deallocate(mpirecvbuf)

call bcast_recv(x%h1x1i,tagh1)
call bcast_recv(x%h2x1i,tagh2)
call bcast_recv(x%h3x1i,tagh3)

call bcast_recv(x%h1x2i,tagh1)      !MZ - these now become x2i sends, and need their own routine in mpimod...
call bcast_recv(x%h2x2i,tagh2)
call bcast_recv(x%h3x2i,tagh3)

call bcast_recv3D_x3i(x%h1x3i,tagh1)
call bcast_recv3D_x3i(x%h2x3i,tagh2)
call bcast_recv3D_x3i(x%h3x3i,tagh3)

call bcast_recv(g1,tagh1)
call bcast_recv(g2,tagh2)
call bcast_recv(g3,tagh3)

call bcast_recv(x%glat,tagglat)
call bcast_recv(x%glon,tagglon)
call bcast_recv(x%alt,tagalt)

call bcast_recv(x%Bmag,tagBmag)
call bcast_recv(x%I,taginc)           !only time that we need to exchange 2D array data, I think
call bcast_recv(x%nullpts,tagnull)    !note that this is not an integer array

allocate(mpirecvbuf(1:lx1,1:lx2,1:lx3))
do icomp=1,3
  call bcast_recv(mpirecvbuf,tageunit1)
  x%e1(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tageunit2)
  x%e2(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tageunit3)
  x%e3(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tager)
  x%er(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tagetheta)
  x%etheta(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tagephi)
  x%ephi(:,:,:,icomp)=mpirecvbuf
end do
deallocate(mpirecvbuf)

call bcast_recv(x%r,tagr)
call bcast_recv(x%theta,tagtheta)
call bcast_recv(x%phi,tagphi)


!COUNT THE NUMBER OF NULL GRID POINTS AND GENERATE A LIST OF NULL INDICES FOR LATER USE
x%lnull=0;
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3)>0.5d0) then
        x%lnull=x%lnull+1
      end if
    end do
  end do
end do
allocate(x%inull(1:x%lnull,1:3))

icount=1
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3)>0.5d0) then
        x%inull(icount,:)=[ix1,ix2,ix3]
        icount=icount+1
      end if
    end do
  end do
end do
print *, 'Done computing null grid points...  Process:  ',myid,' has:  ',x%lnull


!COMPUTE DIFFERENTIAL DISTANCES ALONG EACH DIRECTION
allocate(x%dl1i(1:lx1,1:lx2,1:lx3),x%dl2i(1:lx1,1:lx2,1:lx3),x%dl3i(1:lx1,1:lx2,1:lx3),tmpdx(1:lx1,1:lx2,1:lx3))

tmpdx=spread(spread(x%dx1i,2,lx2),3,lx3)
x%dl1i=tmpdx*x%h1(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx2i,1,lx1),3,lx3)
x%dl2i=tmpdx*x%h2(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx3i,1,lx1),2,lx2)
x%dl3i=tmpdx*x%h3(1:lx1,1:lx2,1:lx3)


deallocate(tmpdx)

end subroutine read_grid_workers


!MZ - need to deallocate x2all-type variables...
subroutine clear_grid(x)

type(curvmesh), intent(inout) :: x

!------------------------------------------------------------
!-------DEALLOCATES GRID VARIABLES.  
!------------------------------------------------------------

deallocate(x%x1,x%x2,x%x3)
deallocate(x%dx1i,x%x1i,x%dx1)
deallocate(x%dx2i,x%x2i,x%dx2)
deallocate(x%dx3i,x%x3i,x%dx3)
deallocate(x%glat,x%glon,x%alt)
deallocate(x%r,x%theta,x%phi)

deallocate(x%h1,x%h2,x%h3)
deallocate(x%h1x1i,x%h2x1i,x%h3x1i)
deallocate(x%h1x2i,x%h2x2i,x%h3x2i)
deallocate(x%h1x3i,x%h2x3i,x%h3x3i)

deallocate(x%I,x%Bmag,x%nullpts)
!    deallocate(x%e1,x%e2,x%e3)    !handled by clear_unitvecs (assuming netural perturbations are used)
!    deallocate(x%er,x%etheta,x%ephi)

deallocate(x%dl1i,x%dl2i,x%dl3i)

if (myid == 0) then
  deallocate(x%x3all,x%x3iall)
  deallocate(x%h1all,x%h2all,x%h3all)
  deallocate(x%h1x1iall,x%h2x1iall,x%h3x1iall)
  deallocate(x%h1x2iall,x%h2x2iall,x%h3x2iall)
  deallocate(x%h1x3iall,x%h2x3iall,x%h3x3iall)
  deallocate(x%rall,x%thetaall,x%phiall)
end if

!THIS NEEDS TO DEALLOCATE BOTH GRAVITY AND THE MAGNETIC FIELD MAGNITUDE
call clear_grav()

end subroutine clear_grid


subroutine clear_unitvecs(x)

!------------------------------------------------------------
!-------DEALLOCATE GRID UNIT VECTORS, WHICH TAKE UP A LOT OF
!-------MEMORY
!------------------------------------------------------------

type(curvmesh), intent(inout) :: x

deallocate(x%e1,x%e2,x%e3,x%er,x%etheta,x%ephi)

end subroutine clear_unitvecs


subroutine load_grav(alt)

!------------------------------------------------------------
!-------LOAD UP GRAV. FIELD ARRAY.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: alt

allocate(g1(lx1,lx2,lx3),g2(lx1,lx2,lx3),g3(lx1,lx2,lx3))

g1=-1*Gconst*Me/(Re+alt)**2

end subroutine load_grav


subroutine clear_grav()
!! DEALLOCATE GRAV. FIELD ARRAY. 

deallocate(g1,g2,g3)

end subroutine clear_grav

end module grid

module grid
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

use mpi, only: mpi_integer, mpi_comm_world, mpi_status_ignore

use phys_consts, only: Gconst,Me,Re,wp,red,black

use mpimod, only: myid, lid, lid2, lid3, &
  tagx1, tagx2, tagx3, tagtheta, tagr, tagphi, tagnull, taglx1, taglx2, taglx3, taglx3all, taginc, &
  tagh1, tagh2, tagh3, tagglat, tagglon, tageunit1, tageunit2, tageunit3, tagetheta, tager, &
  tagalt, tagbmag, tagephi, tagswap, &
  mpi_realprec, taglx2all, tagx3all, tagx2all, &
  bcast_recv, bcast_send, bcast_recv3D_ghost, bcast_send3D_ghost, bcast_recv3D_x3i, bcast_send3D_x3i, &
  bcast_send3D_x2i,bcast_recv3D_x2i, bcast_send1D_2, bcast_recv1D_2, bcast_send1D_3, bcast_recv1D_3

implicit none
private

integer, protected :: lx1,lx2,lx3,lx2all,lx3all    !this is a useful shorthand for most program units using this module, occassionally a program unit needs to define its own size in which case an only statement is required when using this module.

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

  real(wp), dimension(:), allocatable :: x2all
  real(wp), dimension(:), allocatable  :: x2iall
  real(wp), dimension(:), pointer :: dx2all
  real(wp), dimension(:), allocatable  :: dx2iall

  real(wp), dimension(:), allocatable :: x3all
  real(wp), dimension(:), allocatable  :: x3iall
  real(wp), dimension(:), pointer :: dx3all
  real(wp), dimension(:), allocatable  :: dx3iall

  !METRIC FACTORS
  real(wp), dimension(:,:,:), pointer :: h1,h2,h3            !these are the cell-centered metric coefficients
  real(wp), dimension(:,:,:), allocatable :: h1x1i,h2x1i,h3x1i   !metric factors at x1 cell interfaces; dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), allocatable :: h1x2i,h2x2i,h3x2i   !metric factors at x2 interfaces; dim. 2 has lx2+1
  real(wp), dimension(:,:,:), allocatable :: h1x3i,h2x3i,h3x3i   !metric factors at x3 interfaces; dim. 3 has lx3+1

  !ROOT ONLY FULL GRID METRIC FACTORS (WORKERS WILL NOT ALLOCATE)
  real(wp), dimension(:,:,:), pointer :: h1all,h2all,h3all
  real(wp), dimension(:,:,:), allocatable :: h1x1iall,h2x1iall,h3x1iall   !dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), allocatable :: h1x2iall,h2x2iall,h3x2iall   !dim. 2 has lx2all+1
  real(wp), dimension(:,:,:), allocatable :: h1x3iall,h2x3iall,h3x3iall   !dim. 3 has lx3all+1

  !SHALL WE ALSO PRECOMPUTE SOME OF THE PRODUCTS FOR ADVECTION?

  !SIZE INFORMATION
  integer :: lx1,lx2,lx3,lx2all,lx3all    !for program units that may not be able to access module globals

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
  logical :: flagper
end type curvmesh

interface ! read.f90

module subroutine read_grid(indatsize,indatgrid,flagperiodic,x)
character(*), intent(in) :: indatsize,indatgrid
integer, intent(in) :: flagperiodic
type(curvmesh), intent(inout) :: x
end subroutine read_grid

end interface


public :: curvmesh,  lx1,lx2,lx3, lx2all,lx3all, gridflag, flagswap, clear_unitvecs, g1,g2,g3, &
  read_grid, clear_grid, grid_size

contains


subroutine grid_size(indatsize)

!! CHECK THE SIZE OF THE GRID TO BE LOADED AND SET SIZES IN THIS MODULE (NOT IN STRUCTURE THOUGH)

character(*), intent(in) :: indatsize

integer :: iid, ierr, u
logical exists


if (myid==0) then    !root must physically read the size info and pass to workers
  inquire(file=indatsize, exist=exists)
  if (.not.exists) then
     write(stderr,*) 'must generate grid with script before running simulation--grid not present: ',indatsize
     error stop
  endif

  !DETERMINE THE SIZE OF THE GRID TO BE LOADED
  open(newunit=u,file=indatsize,status='old',form='unformatted', &
       access='stream', action='read')
  read(u) lx1,lx2all,lx3all    !note that these are sizes *including ghost cells*
  close(u)
  do iid=1,lid-1
    call mpi_send(lx1,1,MPI_INTEGER,iid,taglx1,MPI_COMM_WORLD,ierr)
    call mpi_send(lx2all,1,MPI_INTEGER,iid,taglx2all,MPI_COMM_WORLD,ierr)
    call mpi_send(lx3all,1,MPI_INTEGER,iid,taglx3all,MPI_COMM_WORLD,ierr)
  end do

  if (ierr/=0) error stop 'grid:grid_size failed mpi_send'

  print *, 'Grid has full size:  ',lx1,lx2all,lx3all
else
  call mpi_recv(lx1,1,MPI_INTEGER,0,taglx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(lx2all,1,MPI_INTEGER,0,taglx2all,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(lx3all,1,MPI_INTEGER,0,taglx3all,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

  if (ierr/=0) error stop 'grid:grid_size failed mpi_recv'
end if

end subroutine grid_size


subroutine clear_grid(x)

type(curvmesh), intent(inout) :: x

!------------------------------------------------------------
!-------DEALLOCATES GRID VARIABLES.
!------------------------------------------------------------

deallocate(x%x3all,x%x2all)
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
  deallocate(x%x2iall,x%dx2all,x%dx2iall)
  deallocate(x%x3iall,x%dx3all,x%dx3iall)
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
!! DEALLOCATE GRID UNIT VECTORS, WHICH TAKE UP A LOT OF MEMORY

type(curvmesh), intent(inout) :: x

deallocate(x%e1,x%e2,x%e3,x%er,x%etheta,x%ephi)

end subroutine clear_unitvecs


subroutine load_grav(alt)
!! LOAD UP GRAV. FIELD ARRAY.  IT IS EXPECTED THAT
!! GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE THEY ARE PASSED INTO THIS ROUTINE.

real(wp), dimension(:,:,:), intent(in) :: alt

allocate(g1(lx1,lx2,lx3),g2(lx1,lx2,lx3),g3(lx1,lx2,lx3))

g1 = -1 * Gconst * Me / (Re + alt)**2

end subroutine load_grav


subroutine clear_grav()
!! DEALLOCATE GRAV. FIELD ARRAY.

deallocate(g1,g2,g3)

end subroutine clear_grav

end module grid

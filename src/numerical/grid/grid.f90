module grid
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

use mesh, only: curvmesh

use phys_consts, only: Gconst,Me,Re,wp,red,black
use reader, only: get_simsize3

use mpimod, only: mpi_integer, mpi_comm_world, mpi_status_ignore, &
myid, lid, lid2, lid3, tag=>gemini_mpi, &
bcast_recv, bcast_send, bcast_recv3D_ghost, bcast_send3D_ghost, bcast_recv3D_x3i, bcast_send3D_x3i, &
bcast_send3D_x2i,bcast_recv3D_x2i, bcast_send1D_2, bcast_recv1D_2, bcast_send1D_3, bcast_recv1D_3

implicit none (type, external)
private
public :: lx1,lx2,lx3, lx2all,lx3all, gridflag, flagswap, clear_unitvecs, g1,g2,g3, &
  read_grid, clear_grid, grid_size, grid_check

external :: mpi_recv, mpi_send

integer, protected :: lx1,lx2,lx3,lx2all,lx3all
!! this is a useful shorthand for most program units using this module,
!! occassionally a program unit needs to define its own size in which case an only statement
!! is required when using this module.

real(wp), dimension(:,:,:), allocatable, protected :: g1,g2,g3
!! gravity, not to be modified by a procedure outside this module
integer, protected :: gridflag
!! for cataloguing the type of grid that we are using, open, closed, inverted, etc.  0 - closed dipole, 1 - inverted open, 2 - standard open.
integer :: flagswap
!! have the x2 and x3 dimensions been swapped?

interface ! read.f90
module subroutine read_grid(indatsize,indatgrid,flagperiodic,x)
character(*), intent(in) :: indatsize,indatgrid
integer, intent(in) :: flagperiodic
type(curvmesh), intent(inout) :: x
end subroutine read_grid
end interface

interface ! check.f90
module subroutine grid_check(x)
type(curvmesh), intent(in) :: x
end subroutine grid_check
end interface

contains


subroutine grid_size(indatsize)
!! CHECK THE SIZE OF THE GRID TO BE LOADED AND SET SIZES IN THIS MODULE (NOT IN STRUCTURE THOUGH)

character(*), intent(in) :: indatsize

if (myid==0) then
  !! root must physically read the size info and pass to workers
  call grid_size_root(indatsize)
else
  call grid_size_worker()
end if

end subroutine grid_size


subroutine grid_size_root(indatsize)
!! DETERMINE THE SIZE OF THE GRID TO BE LOADED
character(*), intent(in) :: indatsize
integer :: i, ierr

call get_simsize3(indatsize, lx1, lx2all, lx3all)

if (lx1 < 1 .or. lx2all < 1 .or. lx3all < 1) then
  write(stderr,*) 'ERROR: reading ' // indatsize
  error stop 'grid_size_root: grid size must be strictly positive'
endif

do i = 1,lid-1
  call mpi_send(lx1,1,MPI_INTEGER, i, tag%lx1,MPI_COMM_WORLD,ierr)
  if (ierr/=0) error stop 'grid_size_root: lx1 failed mpi_send'
  call mpi_send(lx2all,1,MPI_INTEGER, i, tag%lx2all,MPI_COMM_WORLD,ierr)
  if (ierr/=0) error stop 'grid_size_root: lx2all failed mpi_send'
  call mpi_send(lx3all,1,MPI_INTEGER, i, tag%lx3all,MPI_COMM_WORLD,ierr)
  if (ierr/=0) error stop 'grid_size_root: lx3all failed mpi_send'
end do

print *, 'grid_size_root: full grid size:  ',lx1,lx2all,lx3all

end subroutine grid_size_root


subroutine grid_size_worker()

integer :: ierr

call mpi_recv(lx1,1,MPI_INTEGER,0,tag%lx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'grid_size_worker: lx1 failed mpi_send'
call mpi_recv(lx2all,1,MPI_INTEGER,0,tag%lx2all,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'grid_size_worker: lx2all failed mpi_send'
call mpi_recv(lx3all,1,MPI_INTEGER,0,tag%lx3all,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'grid_size_worker: lx3all failed mpi_send'
end subroutine grid_size_worker


subroutine clear_grid(x)
!! DEALLOCATES GRID VARIABLES.

type(curvmesh), intent(inout) :: x

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

if(allocated(x%e1)) deallocate(x%e1,x%e2,x%e3)
if(allocated(x%er)) deallocate(x%er,x%etheta,x%ephi)

deallocate(x%dl1i,x%dl2i,x%dl3i)

if (myid == 0) then
  deallocate(x%x2iall,x%dx2all,x%dx2iall)
  deallocate(x%x3iall,x%dx3all,x%dx3iall)
  deallocate(x%h1all,x%h2all,x%h3all)
  deallocate(x%h1x1iall,x%h2x1iall,x%h3x1iall)
  deallocate(x%h1x2iall,x%h2x2iall,x%h3x2iall)
  deallocate(x%h1x3iall,x%h2x3iall,x%h3x3iall)
  deallocate(x%rall,x%thetaall,x%phiall)
  deallocate(x%altall,x%Bmagall,x%glonall)
end if

!THIS NEEDS TO DEALLOCATE BOTH GRAVITY AND THE MAGNETIC FIELD MAGNITUDE
call clear_grav()

end subroutine clear_grid


subroutine clear_unitvecs(x)
!! DEALLOCATE GRID UNIT VECTORS, WHICH TAKE UP A LOT OF MEMORY

type(curvmesh), intent(inout) :: x

if(allocated(x%e1)) deallocate(x%e1,x%e2,x%e3)
if(allocated(x%er)) deallocate(x%er,x%etheta,x%ephi)

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

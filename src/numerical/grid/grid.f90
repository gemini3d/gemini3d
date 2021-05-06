module grid
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

use meshobj, only: curvmesh
use meshobj_dipole, only: dipolemesh
use meshobj_cart, only: cartmesh
use phys_consts, only: Gconst,Me,Re,wp,red,black
use reader, only: get_simsize3
use mpimod, only: mpi_integer, mpi_comm_world, mpi_status_ignore, &
  mpi_cfg, tag=>gemini_mpi, mpi_realprec, &
  bcast_recv, bcast_send, bcast_recv3D_ghost, bcast_send3D_ghost, bcast_recv3D_x3i, bcast_send3D_x3i, &
  bcast_send3D_x2i,bcast_recv3D_x2i, bcast_send1D_2, bcast_recv1D_2, bcast_send1D_3, bcast_recv1D_3, &
  gather_send3D_ghost,gather_send3D_x2i,gather_send3D_x3i,gather_recv3D_ghost,gather_recv3D_x2i,gather_recv3D_x3i, &
  gather_send,gather_recv,ID2grid

implicit none (type, external)
private
public :: lx1,lx2,lx3, lx2all,lx3all, gridflag, flagswap, g1,g2,g3, &
  read_grid, grid_size, grid_check, grid_drift

external :: mpi_recv, mpi_send

integer, protected :: lx1,lx2,lx3,lx2all,lx3all
!! this is a useful shorthand for most program units using this module,
!! occassionally a program unit needs to define its own size in which case an only statement
!! is required when using this module.
real(wp), dimension(:,:,:), pointer, protected :: g1,g2,g3
!! gravity, not to be modified by a procedure outside this module
integer, protected :: gridflag
!! for cataloguing the type of grid that we are using, open, closed, inverted, etc.  0 - closed dipole, 1 - inverted open, 2 - standard open.
integer :: flagswap
!! have the x2 and x3 dimensions been swapped?

interface ! read.f90
  module subroutine read_grid(indatsize,indatgrid,flagperiodic,x)
    character(*), intent(in) :: indatsize,indatgrid
    integer, intent(in) :: flagperiodic
    class(curvmesh), allocatable, intent(inout) :: x
  end subroutine read_grid
end interface

interface ! check.f90
  module subroutine grid_check(x)
    class(curvmesh), intent(in) :: x
  end subroutine grid_check
end interface

contains


subroutine grid_size(indatsize)
!! CHECK THE SIZE OF THE GRID TO BE LOADED AND SET SIZES IN THIS MODULE (NOT IN STRUCTURE THOUGH)

  character(*), intent(in) :: indatsize
  
  if (mpi_cfg%myid==0) then
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
  
  if (lx1 < 1 .or. lx2all < 1 .or. lx3all < 1) error stop 'grid_size_root: ' // indatsize // ' grid size must be strictly positive'
  
  do i = 1,mpi_cfg%lid-1
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


subroutine grid_drift(x,E02,E03,v2grid,v3grid)
!! Compute the speed the grid is moving at given a background electric field

  class(curvmesh), intent(in) :: x
  reaL(wp), dimension(:,:,:), intent(in) :: E02,E03
  real(wp), intent(out) :: v2grid,v3grid
  integer :: iid,ierr
  real(wp) :: E2ref,E3ref,Bref

  ! Root decides grid drift speed by examining initial background field in this center of its subdomain...
  ! Bad things will happen if these background fields are not uniform, which by definition they should be BUT
  ! no error checking is done on input to insure this, I believe.
  if (mpi_cfg%myid==0) then
    E2ref=E02(lx1,lx2/2,lx3/2)
    E3ref=E03(lx1,lx2/2,lx3/2)
    Bref=x%Bmag(lx1,lx2/2,lx3/2)
    v2grid=E3ref/Bref
    v3grid=-1*E2ref/Bref
    ! FIXME:  error checking to make sure input is sensible for this???
    do iid=1,mpi_cfg%lid-1
      call mpi_send(v2grid,1,mpi_realprec,iid,tag%v2grid,MPI_COMM_WORLD,ierr)
      call mpi_send(v3grid,1,mpi_realprec,iid,tag%v3grid,MPI_COMM_WORLD,ierr)
    end do
  else
    call mpi_recv(v2grid,1,mpi_realprec,0,tag%v2grid,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call mpi_recv(v3grid,1,mpi_realprec,0,tag%v3grid,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  end if
end subroutine grid_drift

end module grid

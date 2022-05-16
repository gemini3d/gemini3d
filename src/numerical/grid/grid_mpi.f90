!> contains procedures for setting up grid that involve message passing of some sort.
module grid_mpi

use, intrinsic:: iso_c_binding, only : C_PTR, c_f_pointer, c_loc, C_INT
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

use meshobj, only: curvmesh
use meshobj_dipole, only: dipolemesh
use meshobj_cart, only: cartmesh
use phys_consts, only: Gconst,Me,Re,wp,red,black
use mpimod, only: mpi_integer, mpi_comm_world, mpi_status_ignore, &
  mpi_cfg, tag=>gemini_mpi, mpi_realprec, &
  bcast_recv, bcast_send, bcast_recv3D_ghost, bcast_send3D_ghost, bcast_recv3D_x3i, bcast_send3D_x3i, &
  bcast_send3D_x2i,bcast_recv3D_x2i, bcast_send1D_2, bcast_recv1D_2, bcast_send1D_3, bcast_recv1D_3, &
  gather_send3D_ghost,gather_send3D_x2i,gather_send3D_x3i,gather_recv3D_ghost,gather_recv3D_x2i,gather_recv3D_x3i, &
  gather_send,gather_recv,ID2grid,grid2ID
use grid, only: lx1,lx2,lx3,lx2all,lx3all,gridflag, &
                get_grid3_coords_raw,get_grid3_coords_hdf5,get_grid3_coords_nc4, &
                set_total_grid_sizes,set_subgrid_sizes,set_gridflag,bind_grav_ptrs

implicit none (type, external)
private
public :: read_grid, grid_check, grid_drift
external :: mpi_recv, mpi_send

interface ! read.f90
  module subroutine read_grid_cart(indatsize,indatgrid,flagperiodic,x, x1,x2,x3,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: indatsize,indatgrid
    integer, intent(in) :: flagperiodic
    class(curvmesh), intent(inout) :: x
    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
    real(wp), intent(in) :: glonctr,glatctr
  end subroutine read_grid_cart

   module subroutine read_grid_dipole(indatsize,indatgrid,flagperiodic,x, x1,x2,x3,x2all,x3all)
    character(*), intent(in) :: indatsize,indatgrid
    integer, intent(in) :: flagperiodic
    class(curvmesh), intent(inout) :: x
    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
  end subroutine read_grid_dipole
end interface

interface ! check.f90
  module subroutine grid_check(x)
    class(curvmesh), intent(in) :: x
  end subroutine grid_check
end interface

contains

!> read in grid and set subgrid sizes; total size must already be set in the grid module via grid_size()
subroutine read_grid(indatsize,indatgrid,flagperiodic, x, xtype, xC)
  character(*), intent(in) :: indatsize,indatgrid
  integer, intent(in) :: flagperiodic
  class(curvmesh), pointer, intent(inout) :: x
  integer(C_INT), intent(inout), optional :: xtype
  type(C_PTR), intent(inout), optional :: xC

  real(wp), dimension(:), allocatable :: x1,x2,x3,x2all,x3all
  integer :: islstart,islfin
  integer, dimension(2) :: indsgrid
  integer iid
  real(wp) :: glonctr,glatctr
  !> For whatever reason c_loc must be called on a static type (not polymorphic) though I don't see why
  !    this limitation exists...
  type(cartmesh), pointer :: xcart
  type(dipolemesh), pointer :: xdipole

  !! Declare grid type that we are dealing with; note lack of matching deallocates assume
  !!   that the compiler will deal with it automatically
  !!  Also set the grid center position if not already dictated by the coordinate system

  call calc_subgrid_size(lx2all,lx3all)
  !! everyone computes what the size of their subgrid should be
  allocate(x1(-1:lx1+2), x2(-1:lx2+2), x3(-1:lx3+2), x2all(-1:lx2all+2), x3all(-1:lx3all+2))
  !! tmp space for coords from file
  call get_grid3_coords(indatgrid,x1,x2all,x3all, glonctr,glatctr)
  !! only need ctr location for certain grid types

  !> each worker needs to set their specific subgrid coordinates
  indsgrid=ID2grid(mpi_cfg%myid, mpi_cfg%lid2)
  !! compute my location on the process grid
  !> x2
  islstart=indsgrid(1)*lx2+1
  !! piece of grid that corresponds to my x3 position
  islfin=islstart+lx2-1
  x2=x2all(islstart-2:islfin+2)
  !> x3
  islstart=indsgrid(2)*lx3+1
  !! piece of grid that corresponds to my x3 position
  islfin=islstart+lx3-1
  x3=x3all(islstart-2:islfin+2)


  print*, 'read_grid has size:  ',lx1,lx2,lx3,lx2all,lx3all
  !! FIXME: hardcode grid type for now; compute it from the coordinates eventually??
  !! right now we just have Cartesian and dipole so it's easy to detect based on x2
  if (maxval(abs(x2))<100) then
    print*, ' Detected dipole grid...'
    !allocate(dipolemesh::x)
    allocate(xdipole)
    x=>xdipole
    call read_grid_dipole(indatsize,indatgrid,flagperiodic,x,x1,x2,x3,x2all,x3all)
    if (present(xC) .and. present(xtype)) then
      xC = c_loc(xdipole)
      xtype=2
      print *, "dipole c_loc"
    end if
  else
    print*, 'Detected Cartesian grid...'
    !allocate(cartmesh::x)
    allocate(xcart)
    x=>xcart
    call read_grid_cart(indatsize,indatgrid,flagperiodic,x,x1,x2,x3,x2all,x3all,glonctr,glatctr)
    print*, 'read_grid_cart done'
    if (present(xC) .and. present(xtype)) then
      xC = c_loc(xcart)
      xtype=1
      print *, "xcart c_loc"
    end if
  end if
  print*, 'read_grid end has size:  ',lx1,lx2,lx3,lx2all,lx3all
  print*, 'grid object has size:  ',x%lx1,x%lx2,x%lx3
end subroutine read_grid

subroutine get_grid3_coords(path,x1,x2all,x3all,glonctr,glatctr)
  character(*), intent(in) :: path
  real(wp), dimension(:), intent(inout) :: x1,x2all,x3all
  real(wp) :: glonctr,glatctr

  character(:), allocatable :: fmt

  fmt = path(index(path, '.', back=.true.) : len(path))
  select case (fmt)
    case ('.dat')
      call get_grid3_coords_raw(path,x1,x2all,x3all,glonctr,glatctr)
    case ('.h5')
      call get_grid3_coords_hdf5(path,x1,x2all,x3all,glonctr,glatctr)
    case ('.nc')
      call get_grid3_coords_nc4(path,x1,x2all,x3all,glonctr,glatctr)
    case default
      error stop 'grid:read:get_grid3: unknown grid format: ' // fmt
  end select

  if(size(x1) < 1) error stop 'grid:get_grid3_coords: size(x1) must be strictly positive'
  if(size(x2all) < 1) error stop 'grid:get_grid3_coords: size(x2all) must be strictly positive'
  if(size(x3all) < 1) error stop 'grid:get_grid3_coords: size(x3all) must be strictly positive'
end subroutine get_grid3_coords


!> worker subgrid sizes; requires knowledge of mpi, though not any direct mpi calls
subroutine calc_subgrid_size(lx2all, lx3all)
  integer, intent(in) :: lx2all, lx3all
  integer :: lx2, lx3

  !! use only non-swapped axes
  if(lx2all==1) then
    print *, 'get_subgrid_size: 2D run with singleton x2'
    lx2 = 1
    lx3 = lx3all/mpi_cfg%lid
  else if (lx3all==1) then
    print*, 'get_subgrid_size:  2D run with singleton x3'
    lx3=1
    lx2=lx2all/mpi_cfg%lid
  else
    print *, 'get_subgrid_size: 3D run'
    !! should divide evenly if generated from process_grid
    lx2 = lx2all/mpi_cfg%lid2
    lx3 = lx3all/mpi_cfg%lid3
  end if

  if(lx1 < 1) error stop 'grid:calc_subgrid_size: lx1 must be strictly positive'
  if(lx2 < 1) error stop 'grid:calc_subgrid_size: lx2 must be strictly positive'
  if(lx3 < 1) error stop 'grid:calc_subgrid_size: lx3 must be strictly positive'
  if(lx2all < lx2) error stop 'grid:calc_subgrid_size: lx2all must be > lx2'
  if(lx3all < lx3) error stop 'grid:calc_subgrid_size: lx3all must be > lx3'

  if(lx2all > 1 .and. lx3all > 1) then
    if(lx2 == 1 .or. lx3 == 1) error stop "read_grid_root: 3D grids cannot be partitioned with a single MPI image on an axis"
  end if

  call set_subgrid_sizes(lx2,lx3)
end subroutine calc_subgrid_size


!> compute grid drift speed
subroutine grid_drift(x,E02,E03,v2grid,v3grid)
!! Compute the speed the grid is moving at given a background electric field
  class(curvmesh), intent(in) :: x
  reaL(wp), dimension(:,:,:), intent(in) :: E02,E03
  real(wp), intent(inout) :: v2grid,v3grid
  !! intent(out)
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

end module grid_mpi

!> contains procedures for setting up grid that involve message passing of some sort.
submodule (grid) grid_mpi

use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

use phys_consts, only: Gconst,Me,Re,red,black
use mpimod, only: mpi_cfg, tag=>gemini_mpi, mpi_realprec, &
  bcast_recv, bcast_send, bcast_recv3D_ghost, bcast_send3D_ghost, bcast_recv3D_x3i, bcast_send3D_x3i, &
  bcast_send3D_x2i,bcast_recv3D_x2i, bcast_send1D_2, bcast_recv1D_2, bcast_send1D_3, bcast_recv1D_3, &
  gather_send3D_ghost,gather_send3D_x2i,gather_send3D_x3i,gather_recv3D_ghost,gather_recv3D_x2i,gather_recv3D_x3i, &
  gather_send,gather_recv,ID2grid,grid2ID

use mpi_f08, only: mpi_integer, mpi_comm_world, mpi_status_ignore,mpi_recv, mpi_send

implicit none (type, external)

interface ! read.f90
  module subroutine read_grid_cartdip(indatsize,indatgrid,flagperiodic,x, x1,x2,x3,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: indatsize,indatgrid
    integer, intent(in) :: flagperiodic
    class(curvmesh), intent(inout) :: x
    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
    real(wp), intent(in) :: glonctr,glatctr
  end subroutine read_grid_cartdip
!  module subroutine read_grid_dipole(indatsize,indatgrid,flagperiodic,x, x1,x2,x3,x2all,x3all)
!    character(*), intent(in) :: indatsize,indatgrid
!    integer, intent(in) :: flagperiodic
!    class(curvmesh), intent(inout) :: x
!    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
!  end subroutine read_grid_dipole
end interface

contains

  module procedure read_grid

    real(wp), dimension(:), allocatable :: x2,x3,x2all,x3all
    integer :: islstart,islfin
    integer, dimension(2) :: indsgrid
    real(wp) :: glonctr,glatctr
    !> For whatever reason c_loc must be called on a static type (not polymorphic) though I don't see why
    !    this limitation exists...
    ! type(cartmesh), pointer :: xcart
    ! type(dipolemesh), pointer :: xdipole
    ! integer :: gridtype

    !! Declare grid type that we are dealing with; note lack of matching deallocates assume
    !!   that the compiler will deal with it automatically
    !!  Also set the grid center position if not already dictated by the coordinate system

    !call calc_subgrid_size(lx2all,lx3all)
    !! everyone computes what the size of their subgrid should be
    !^ this is now done as a separate step from the main application through libgemini_mpi calls
    call alloc_x1coords(lx1)
    allocate(x2(-1:lx2+2), x3(-1:lx3+2), x2all(-1:lx2all+2), x3all(-1:lx3all+2))
    !! tmp space for coords from file
    call get_grid3_coords(indatgrid,x1,x2all,x3all, glonctr,glatctr)
    call set_fullgrid_lims(x2all,x3all)
    !! read the grid coordinates in from a file only need ctr location for certain grid types

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

    !> allocate the correct class for the grid and bind a C pointer (which can only be done at creation!)
    if (present(xC) .and. present(xtype)) then
      call meshobj_alloc(x1,x2,x3,x,xtype,xC)
    else
      call meshobj_alloc(x1,x2,x3,x)
    end if

    !> execute call to collect global data, as needed
    call read_grid_cartdip(indatsize,indatgrid,flagperiodic,x,x1,x2,x3,x2all,x3all,glonctr,glatctr)    ! dipole grid doesn't use ctr coords

    !> just to be careful explicitly deallocate temp arrays
    deallocate(x2,x3,x2all,x3all)
  end procedure read_grid


  module procedure calc_subgrid_size

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
  end procedure calc_subgrid_size


  module procedure grid_drift
    integer :: iid
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
        call mpi_send(v2grid,1,mpi_realprec,iid,tag%v2grid,MPI_COMM_WORLD)
        call mpi_send(v3grid,1,mpi_realprec,iid,tag%v3grid,MPI_COMM_WORLD)
      end do
    else
      call mpi_recv(v2grid,1,mpi_realprec,0,tag%v2grid,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      call mpi_recv(v3grid,1,mpi_realprec,0,tag%v3grid,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
    end if
  end procedure grid_drift

end submodule grid_mpi

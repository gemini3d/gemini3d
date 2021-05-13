submodule (grid) grid_read

!use mpimod, only : mpi_realprec

implicit none (type, external)

interface ! readgrid_*.f90
  module subroutine get_grid3_coords_raw(path,x1,x2all,x3all)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(out) :: x1,x2all,x3all
  end subroutine get_grid3_coords_raw

  module subroutine get_grid3_coords_hdf5(path,x1,x2all,x3all)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(out) :: x1,x2all,x3all
  end subroutine get_grid3_coords_hdf5

  module subroutine get_grid3_coords_nc4(path,x1,x2all,x3all)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(out) :: x1,x2all,x3all
  end subroutine get_grid3_coords_nc4
end interface

contains

!> Read in the grid information and prep grid object.  Note that there are also module-scope variables
!   that are (redundantly, for convenience) defined based on the grid structure and this procedure 
!   must also set those.  
module procedure read_grid
! subroutine read_grid(indatsize,indatgrid,flagperiodic,x)
  real(wp), dimension(:), allocatable :: x1,x2,x3,x2all,x3all
  integer :: gridtype
  integer :: islstart,islfin
  integer, dimension(2) :: indsgrid
  integer iid

  call set_subgrid_size()    ! everyone computes what the size of their subgrid should be
  allocate(x1(-1:lx1+2),x2(-1:lx2+2),x3(-1:lx3+2),x2all(-1:lx2all+2),x3all(-1:lx3all+2))   ! tmp space for coords from file
  call get_grid3_coords(indatgrid,x1,x2all,x3all)

  !> each worker needs to set their specific subgrid coordinates
  indsgrid=ID2grid(mpi_cfg%myid, mpi_cfg%lid2)     !compute my location on the process grid
  !! x2
  islstart=indsgrid(1)*lx2+1              !piece of grid that corresponds to my x3 position
  islfin=islstart+lx2-1
  x2=x2all(islstart-2:islfin+2)
  !! x3
  islstart=indsgrid(2)*lx3+1              !piece of grid that corresponds to my x3 position
  islfin=islstart+lx3-1
  x3=x3all(islstart-2:islfin+2)

  ! FIXME: hardcode grid type for now; compute it from the coordinates eventually??
  gridtype=1
  
  !> Declare grid type that we are dealing with; note lack of matching deallocates assume
  !   that the compiler will deal with it automatically
  select case (gridtype)
    case(0)    ! cartesian
      allocate(cartmesh::x) 
    case(1)    ! dipole
      allocate(dipolemesh::x)
    case default
      error stop ' invalid mesh type specified; cannot instantiate!'
  end select

  !> Create the grid object
  call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays
  deallocate(x1,x2,x3,x2all,x3all)
  call x%init()                              ! allocate space for subgrid variables
  call x%make()                              ! fill auxiliary arrays
  
  !> We need to collect the info for root's fullgrid variables
  print*, 'subgrid metric',mpi_cfg%myid,minval(x%h2x1i),maxval(x%h2x1i)
  print*, 'subgrid theta',mpi_cfg%myid,minval(x%theta),maxval(x%theta)

  print*, '...collecting full mesh variables...'
  if (mpi_cfg%myid==0) then
    call x%init_storage_root()                ! now we have space in type to store full-grid arrays for gather
    call gather_grid_root(x%h1,x%h2,x%h3, &
                          x%h1x1i,x%h2x1i,x%h3x1i, &
                          x%h1x2i,x%h2x2i,x%h3x2i, &
                          x%h1x3i,x%h2x3i,x%h3x3i, &
                          x%r,x%theta,x%phi, &
                          x%alt,x%Bmag,x%glon, &
                          x%h1all,x%h2all,x%h3all, &
                          x%h1x1iall,x%h2x1iall,x%h3x1iall, &
                          x%h1x2iall,x%h2x2iall,x%h3x2iall, &
                          x%h1x3iall,x%h2x3iall,x%h3x3iall, &
                          x%rall,x%thetaall,x%phiall, &
                          x%altall,x%Bmagall,x%glonall)     
    !! note that we can fill arrays manually with our own routines rather than use x%set_root, saves temp arrays and memory
    call x%calc_coord_diffs_root()
    !print*, 'Fullgrid tests...'
    !print*, minval(x%h1all),maxval(x%h1all), shape(x%h1all)
    !print*, minval(x%h2all),maxval(x%h2all), shape(x%h2all)
    !print*, minval(x%h3all),maxval(x%h3all), shape(x%h3all)

    print*, minval(x%h1x1iall),maxval(x%h1x1iall), shape(x%h1x1iall)
    print*, minval(x%h2x1iall),maxval(x%h2x1iall), shape(x%h2x1iall)
    print*, minval(x%h3x1iall),maxval(x%h3x1iall), shape(x%h3x1iall)

    !print*, minval(x%h1x2iall),maxval(x%h1x2iall)
    !print*, minval(x%h2x2iall),maxval(x%h2x2iall)
    !print*, minval(x%h3x2iall),maxval(x%h3x2iall)

    !print*, minval(x%h1x3iall),maxval(x%h1x3iall)
    !print*, minval(x%h2x3iall),maxval(x%h2x3iall)
    !print*, minval(x%h3x3iall),maxval(x%h3x3iall)

    !print*, minval(x%rall),maxval(x%rall)
    !print*, minval(x%thetaall),maxval(x%thetaall)
    !print*, minval(x%phiall),maxval(x%phiall)

    !print*, minval(x%altall),maxval(x%altall)
    !print*, minval(x%Bmagall),maxval(x%Bmagall)
    !print*, minval(x%glonall),maxval(x%glonall)
  else
    !! gather
    call gather_grid_workers(x%h1,x%h2,x%h3, &
                          x%h1x1i,x%h2x1i,x%h3x1i, &
                          x%h1x2i,x%h2x2i,x%h3x2i, &
                          x%h1x3i,x%h2x3i,x%h3x3i, &
                          x%r,x%theta,x%phi, &
                          x%alt,x%Bmag,x%glon)
  end if
  
  !> Assign periodic or not based on user input
  call x%set_periodic(flagperiodic)
  
  !> Set flags for module scope vars.
  gridflag=x%gridflag
  
  !> Set gravitational fields for module scope vars., use pointers to avoid duplicating data
  g1=>x%g1; g2=>x%g2; g3=>x%g3
  
  !> Make sure we have a sensible x2,3 decomposition of grid
  !> and that parameters aren't impossible
  if(mpi_cfg%myid == 0) call grid_check(x)
end procedure read_grid


subroutine set_subgrid_size()

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

  ! FIXME: right now just force this to zero so later swap-specific code does not get triggered (eventually needs to be removed)
  flagswap=0

  if(lx2all > 1 .and. lx3all > 1) then
    if(lx2 == 1 .or. lx3 == 1) error stop "read_grid_root: 3D grids cannot be partitioned with a single MPI image on an axis"
  end if
end subroutine set_subgrid_size


subroutine get_grid3_coords(path,x1,x2all,x3all)
  character(*), intent(in) :: path
  real(wp), dimension(:), intent(out) :: x1,x2all,x3all

  character(:), allocatable :: fmt

  fmt = path(index(path, '.', back=.true.) : len(path))
  select case (fmt)
    case ('.dat')
      call get_grid3_coords_raw(path,x1,x2all,x3all)
    case ('.h5')
      call get_grid3_coords_hdf5(path,x1,x2all,x3all)
    case ('.nc')
      call get_grid3_coords_nc4(path,x1,x2all,x3all)
    case default
      write(stderr,*) 'grid:read:get_grid3: unknown grid format: ' // fmt
      error stop 2
  end select 
end subroutine get_grid3_coords


!> pull full grid vars. from workers into root arrays
subroutine gather_grid_root(h1,h2,h3, &
                        h1x1i,h2x1i,h3x1i, &
                        h1x2i,h2x2i,h3x2i, &
                        h1x3i,h2x3i,h3x3i, &
                        r,theta,phi, &
                        alt,Bmag,glon, &
                        h1all,h2all,h3all, &
                        h1x1iall,h2x1iall,h3x1iall, &
                        h1x2iall,h2x2iall,h3x2iall, &
                        h1x3iall,h2x3iall,h3x3iall, &
                        rall,thetaall,phiall, &
                        altall,Bmagall,glonall)
  real(wp), dimension(:,:,:), intent(in) :: h1,h2,h3
  real(wp), dimension(:,:,:), intent(in) :: h1x1i,h2x1i,h3x1i
  real(wp), dimension(:,:,:), intent(in) :: h1x2i,h2x2i,h3x2i
  real(wp), dimension(:,:,:), intent(in) :: h1x3i,h2x3i,h3x3i
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(:,:,:), intent(in) :: alt,Bmag,glon
  real(wp), dimension(:,:,:), intent(out) :: h1all,h2all,h3all
  real(wp), dimension(:,:,:), intent(out) :: h1x1iall,h2x1iall,h3x1iall
  real(wp), dimension(:,:,:), intent(out) :: h1x2iall,h2x2iall,h3x2iall
  real(wp), dimension(:,:,:), intent(out) :: h1x3iall,h2x3iall,h3x3iall
  real(wp), dimension(:,:,:), intent(out) :: rall,thetaall,phiall
  real(wp), dimension(:,:,:), intent(out) :: altall,Bmagall,glonall

  call gather_recv3D_ghost(h1,tag%h1,h1all)
  call gather_recv3D_ghost(h2,tag%h2,h2all)
  call gather_recv3D_ghost(h3,tag%h3,h3all)

  print*, 1

  call gather_recv(h1x1i,tag%h1,h1x1iall)
  call gather_recv(h2x1i,tag%h2,h2x1iall)
  call gather_recv(h3x1i,tag%h3,h3x1iall)

  print*, 2

  call gather_recv3D_x2i(h1x2i,tag%h1,h1x2iall)
  call gather_recv3D_x2i(h2x2i,tag%h2,h2x2iall)
  call gather_recv3D_x2i(h3x2i,tag%h3,h3x2iall)

  print*, 3

  call gather_recv3D_x3i(h1x3i,tag%h1,h1x3iall)
  call gather_recv3D_x3i(h2x3i,tag%h2,h2x3iall)
  call gather_recv3D_x3i(h3x3i,tag%h3,h3x3iall)

  print*, 4

  call gather_recv(r,tag%r,rall)
  call gather_recv(theta,tag%theta,thetaall)
  call gather_recv(phi,tag%phi,phiall)

  print*, 5

  call gather_recv(alt,tag%alt,altall)
  call gather_recv(Bmag,tag%Bmag,Bmagall)
  call gather_recv(glon,tag%glon,glonall)

  print*, 6
end subroutine gather_grid_root


!> send full grid vars. to root
subroutine gather_grid_workers(h1,h2,h3, &
                        h1x1i,h2x1i,h3x1i, &
                        h1x2i,h2x2i,h3x2i, &
                        h1x3i,h2x3i,h3x3i, &
                        r,theta,phi, &
                        alt,Bmag,glon)
  real(wp), dimension(:,:,:), intent(in) :: h1,h2,h3
  real(wp), dimension(:,:,:), intent(in) :: h1x1i,h2x1i,h3x1i
  real(wp), dimension(:,:,:), intent(in) :: h1x2i,h2x2i,h3x2i
  real(wp), dimension(:,:,:), intent(in) :: h1x3i,h2x3i,h3x3i
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(:,:,:), intent(in) :: alt,Bmag,glon

  call gather_send3D_ghost(h1,tag%h1)
  call gather_send3D_ghost(h2,tag%h2)
  call gather_send3D_ghost(h3,tag%h3)

  call gather_send(h1x1i,tag%h1)
  call gather_send(h2x1i,tag%h2)
  call gather_send(h3x1i,tag%h3)

  call gather_send3D_x2i(h1x2i,tag%h1)
  call gather_send3D_x2i(h2x2i,tag%h2)
  call gather_send3D_x2i(h3x2i,tag%h3)

  call gather_send3D_x3i(h1x3i,tag%h1)
  call gather_send3D_x3i(h2x3i,tag%h2)
  call gather_send3D_x3i(h3x3i,tag%h3)

  call gather_send(r,tag%r)
  call gather_send(theta,tag%theta)
  call gather_send(phi,tag%phi)

  call gather_send(alt,tag%alt)
  call gather_send(Bmag,tag%Bmag)
  call gather_send(glon,tag%glon)
end subroutine gather_grid_workers

end submodule grid_read

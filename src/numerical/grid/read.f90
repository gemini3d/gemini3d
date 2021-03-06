submodule (grid) grid_read

!use mpimod, only : mpi_realprec

implicit none (type, external)

contains

!> Read in the grid information and prep grid object.  Note that there are also module-scope variables
!   that are (redundantly, for convenience) defined based on the grid structure and this procedure
!   must also set those.
module procedure read_grid_cart
! subroutine read_grid(indatsize,indatgrid,flagperiodic,x)

  call x%set_center(glonctr,glatctr)

  !> Create the grid object
  call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays

  call x%init()                              ! allocate space for subgrid variables
  call x%make()                              ! fill auxiliary arrays

  !> We need to collect the info for root's fullgrid variables
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
end procedure read_grid_cart


module procedure read_grid_dipole
! subroutine read_grid(indatsize,indatgrid,flagperiodic,x)

  !> Create the grid object
  call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays

  call x%init()                              ! allocate space for subgrid variables
  call x%make()                              ! fill auxiliary arrays

  !> We need to collect the info for root's fullgrid variables
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
end procedure read_grid_dipole

!--------------------------------------------------------------------------------------------------



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
  real(wp), dimension(:,:,:), intent(inout) :: h1all,h2all,h3all
  real(wp), dimension(:,:,:), intent(inout) :: h1x1iall,h2x1iall,h3x1iall
  real(wp), dimension(:,:,:), intent(inout) :: h1x2iall,h2x2iall,h3x2iall
  real(wp), dimension(:,:,:), intent(inout) :: h1x3iall,h2x3iall,h3x3iall
  real(wp), dimension(:,:,:), intent(inout) :: rall,thetaall,phiall
  real(wp), dimension(:,:,:), intent(inout) :: altall,Bmagall,glonall

  call gather_recv3D_ghost(h1,tag%h1,h1all)
  call gather_recv3D_ghost(h2,tag%h2,h2all)
  call gather_recv3D_ghost(h3,tag%h3,h3all)

  call gather_recv(h1x1i,tag%h1,h1x1iall)
  call gather_recv(h2x1i,tag%h2,h2x1iall)
  call gather_recv(h3x1i,tag%h3,h3x1iall)

  call gather_recv3D_x2i(h1x2i,tag%h1,h1x2iall)
  call gather_recv3D_x2i(h2x2i,tag%h2,h2x2iall)
  call gather_recv3D_x2i(h3x2i,tag%h3,h3x2iall)

  call gather_recv3D_x3i(h1x3i,tag%h1,h1x3iall)
  call gather_recv3D_x3i(h2x3i,tag%h2,h2x3iall)
  call gather_recv3D_x3i(h3x3i,tag%h3,h3x3iall)

  call gather_recv(r,tag%r,rall)
  call gather_recv(theta,tag%theta,thetaall)
  call gather_recv(phi,tag%phi,phiall)

  call gather_recv(alt,tag%alt,altall)
  call gather_recv(Bmag,tag%Bmag,Bmagall)
  call gather_recv(glon,tag%glon,glonall)

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

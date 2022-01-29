!> drift-related calculations for libgemini_mpi
submodule (gemini3d_mpi) libgem_mpi_drifts

implicit none (type, external)

contains
  !> add in background field, accounting for whether the user specified a lagrangian grid
!  module subroutine BGfield_Lagrangian(cfg,x,v2grid,v3grid,E1,E2,E3)
!    type(gemini_cfg), intent(in) :: cfg
!    class(curvmesh), intent(in) :: x
!    real(wp), intent(inout) :: v2grid,v3grid
!    real(wp), dimension(:,:,:), intent(inout) :: E1,E2,E3
!  end subroutine BGfield_Lagrangian
  module procedure BGfield_Lagrangian
    real(wp), dimension(:,:,:), allocatable :: E01,E02,E03
    integer :: lx1,lx2,lx3

    lx1=size(E2,1); lx2=size(E2,2); lx3=size(E3,3);
    allocate(E01(lx1,lx2,lx3),E02(lx1,lx2,lx3),E03(lx1,lx2,lx3))
    E01=0; E02=0; E03=0;
    if (cfg%flagE0file==1) then
      call get_BGEfields(x,E01,E02,E03)
    end if
    if (cfg%flaglagrangian) then    ! Lagrangian (moving) grid; compute from input background electric fields
      call grid_drift(x,E02,E03,v2grid,v3grid)
      if (mpi_cfg%myid==0) print*, mpi_cfg%myid,' using Lagrangian grid moving at:  ',v2grid,v3grid
    else                            ! stationary grid
      v2grid = 0
      v3grid = 0
      E1 = E1 + E01    ! FIXME: this is before dist fields are computed???
      E2 = E2 + E02
      E3 = E3 + E03
    end if
    deallocate(E01,E02,E03)
  end procedure BGfield_Lagrangian


  !> Compute initial perp drifts
!  module subroutine get_initial_drifts(cfg,x,nn,Tn,vn1,vn2,vn3,ns,Ts,vs1,vs2,vs3,B1,E2,E3)
!    type(gemini_cfg), intent(in) :: cfg
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:,:), intent(in) :: nn
!    real(wp), dimension(:,:,:), intent(in) :: Tn,vn1,vn2,vn3
!    real(wp), dimension(:,:,:,:), intent(in) :: ns,Ts,vs1
!    real(wp), dimension(:,:,:,:), intent(inout) :: vs2,vs3
!    real(wp), dimension(:,:,:), intent(in) :: B1
!    real(wp), dimension(:,:,:), intent(in) :: E2,E3
!  end subroutine get_initial_drifts
  module procedure get_initial_drifts
    real(wp), dimension(:,:,:), allocatable :: sig0,sigP,sigH,sigPgrav,sigHgrav
    real(wp), dimension(:,:,:,:), allocatable :: muP,muH,nusn
    integer :: lx1,lx2,lx3,lsp
    
    lx1=x%lx1; lx2=x%lx2; lx3=x%lx3; lsp=size(ns,4);
    allocate(sig0(lx1,lx2,lx3),sigP(lx1,lx2,lx3),sigH(lx1,lx2,lx3),sigPgrav(lx1,lx2,lx3),sigHgrav(lx1,lx2,lx3))
    allocate(muP(lx1,lx2,lx3,lsp),muH(lx1,lx2,lx3,lsp),nusn(lx1,lx2,lx3,lsp))
    call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
    call velocities(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)
    deallocate(sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
  end procedure get_initial_drifts
end submodule libgem_mpi_drifts

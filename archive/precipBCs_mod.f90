!
!  subroutine precipBCs_GDIKHI(t,x,W0,PhiWmWm2)
!
!    !------------------------------------------------------------
!    !-------LOAD UP ARRAYS CONTAINING TOP BOUNDARY CHAR. ENERGY
!    !-------AND TOTAL ENERGY FLUX.  GRID VARIABLES INCLUDE
!    !-------GHOST CELLS
!    !------------------------------------------------------------
!
!    real(wp), intent(in) :: t
!    type(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:), intent(out) :: W0,PhiWmWm2
!
!    real(wp) :: W0pk,PhiWpk,meanW0x3,meanPhiWx3,sigW0x3,sigPhiWx3
!    real(wp) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve,sigt,meant
!    integer :: ix2,ix3,iprec,lx2,lx3,lprec
!
!    
!    lx2=size(W0,1)
!    lx3=size(W0,2)
!    lprec=size(W0,3)    !assumed to be 2 in this subroutine
!
!
!    !BACKGROUND PRECIPITATION
!    W0pk=3d3
!    PhiWpk=1d-5
!    do ix3=1,lx3
!      do ix2=1,lx2
!        W0(ix2,ix3,1)=W0pk 
!        PhiWmWm2(ix2,ix3,1)=PhiWpk       
!      end do
!    end do
!
!
!    !PARAMETERS FOR DISTURBANCE PRECIPITATION
!    W0pk=100d0
!!      sigW0x3=100d3
!!      meanW0x3=0d0
!    PhiWpk=0d0
!!      PhiWpk=1d-5    !successful grad-drift attempts
!!      PhiWpk=1d-4    !Swoboda blur testing
!!      PhiWpk=0.05d0    !testing of convergent Hall drifts
!!      PhiWpk=5d0
!!      sigPhiWx3=100d3
!!      meanPhiWx3=0d0
!
!!      W0pk=0.3d3
!!      sigW0x3=100d3
!!      meanW0x3=0d0
!!      PhiWpk=2d0
!!      sigPhiWx3=100d3
!!      meanPhiWx3=0d0
!
!    sigx2=50d3
!    meanx2=0d0
!!    sigx3=10d3
!    sigx3=25d3
!    meant=900d0
!    sigt=450d0
!    x30amp=0d3
!    varc=200d0
!
!    !DISTURBANCE ELECTRON PRECIPITATION PATTERN
!    do ix3=1,lx3
!      do ix2=1,lx2
!        W0(ix2,ix3,2)=W0pk
!        PhiWmWm2(ix2,ix3,2)=PhiWpk
!      end do
!    end do
!
!  end subroutine precipBCs_GDIKHI
!
!
!  subroutine precipBCs_3DPCarc2(t,x,W0,PhiWmWm2)
!
!    !------------------------------------------------------------
!    !-------LOAD UP ARRAYS CONTAINING TOP BOUNDARY CHAR. ENERGY
!    !-------AND TOTAL ENERGY FLUX.  GRID VARIABLES INCLUDE
!    !-------GHOST CELLS
!    !------------------------------------------------------------
!
!    real(wp), intent(in) :: t
!    type(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:), intent(out) :: W0,PhiWmWm2
!
!    real(wp) :: W0pk,PhiWpk,meanW0x3,meanPhiWx3,sigW0x3,sigPhiWx3
!    real(wp) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve,sigt,meant
!    integer :: ix2,ix3,iprec,lprec
!
!    
!    lprec=size(W0,3)    !assumed to be 2 in this subroutine
!
!
!    !BACKGROUND PRECIPITATION
!    W0pk=5000d0
!    PhiWpk=0.01d0
!    do ix3=1,lx3
!      do ix2=1,lx2
!        W0(ix2,ix3,1)=W0pk
!        PhiWmWm2(ix2,ix3,1)=PhiWpk
!      end do
!    end do
!
!
!    !PARAMETERS FOR DISTURBANCE PRECIPITATION
!    W0pk=100d0
!!      sigW0x3=100d3
!!      meanW0x3=0d0
!    PhiWpk=0.4d0
!    sigx2=50d3
!    meanx2=0d0
!!    sigx3=10d3
!    sigx3=25d3
!    meant=900d0
!    sigt=450d0
!    x30amp=0d3
!    varc=200d0
!
!    !DISTURBANCE ELECTRON PRECIPITATION PATTERN
!    do ix3=1,lx3
!      do ix2=1,lx2
!        meanx3=-varc*meant+varc*t
!
!        W0(ix2,ix3,2)=W0pk
!        PhiWmWm2(ix2,ix3,2)=(PhiWpk*exp(-(x%x3(ix3)-(-sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x%x2(ix2)) ))**2/2d0/sigx3**2)- &
!                          Phiwpk*exp(-(x%x3(ix3)-(sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x%x2(ix2)) ))**2/2d0/sigx3**2) )* &
!                            exp(-(x%x2(ix2)-meanx2)**2/2d0/sigx2**2)*exp(-(t-meant)**2/2d0/sigt**2)
!        if (PhiWmWm2(ix2,ix3,2) < 0d0) then
!          PhiWmWm2(ix2,ix3,2)=0d0
!        end if
!      end do
!    end do
!
!  end subroutine precipBCs_3DPCarc2
!
!
!  subroutine precipBCs_original(t,x,W0,PhiWmWm2)
!
!    !------------------------------------------------------------
!    !-------LOAD UP ARRAYS CONTAINING TOP BOUNDARY CHAR. ENERGY
!    !-------AND TOTAL ENERGY FLUX.  GRID VARIABLES INCLUDE
!    !-------GHOST CELLS
!    !------------------------------------------------------------
!
!    real(wp), intent(in) :: t
!    type(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:), intent(out) :: W0,PhiWmWm2
!
!    real(wp) :: W0pk,PhiWpk,meanW0x3,meanPhiWx3,sigW0x3,sigPhiWx3
!    real(wp) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve,sigt,meant
!    integer :: ix2,ix3,iprec,lprec
!
!    
!!    lx2=size(W0,1)
!!    lx3=size(W0,2)
!    lprec=size(W0,3)    !assumed to be 2 in this subroutine (background and distubance precipitation)
!
!
!    !BACKGROUND PRECIPITATION
!    W0pk=3d3
!    PhiWpk=0.1d0
!    do ix3=1,lx3
!      do ix2=1,lx2
!        W0(ix2,ix3,1)=W0pk 
!        PhiWmWm2(ix2,ix3,1)=PhiWpk       
!      end do
!    end do
!
!
!    !PARAMETERS FOR DISTURBANCE PRECIPITATION
!    W0pk=300d0
!!    if (t<1800) then
!!      PhiWpk=1d-5    !something very small, but nonzero
!!    else
!!      PhiWpk=2d0
!!    end if
!    PhiWpk=1d-5
!    sigx2=0.33d0*(x%x2(lx2)-x%x2(1))
!    meanx2=0.5d0*(x%x2(1)+x%x2(lx2))
!!    sigx3=0.33d0*(x3(lx3)-x3(1))     !needs to be calculated based on full grid
!!    meanx3=0.5d0*(x3(1)+x3(lx3))     !ditto
!    sigx3=0.33d0*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
!    meanx3=0.5d0*(x%x3all(1)+x%x3all(lx3all))
!
!
!    !DISTURBANCE ELECTRON PRECIPITATION PATTERN
!    do ix3=1,lx3
!      do ix2=1,lx2
!        W0(ix2,ix3,2)=W0pk
!        PhiWmWm2(ix2,ix3,2)=PhiWpk*exp(-(x%x3(ix3)-meanx3)**2/2d0/sigx3**2)* &
!                            exp(-(x%x2(ix2)-meanx2)**2/2d0/sigx2**2)
!        if (PhiWmWm2(ix2,ix3,2) < 0d0) then
!          PhiWmWm2(ix2,ix3,2)=0d0
!        end if
!      end do
!    end do
!
!
!!    !PARAMETERS FOR DISTURBANCE PRECIPITATION
!!    W0pk=100d0
!!!      sigW0x3=100d3
!!!      meanW0x3=0d0
!!!    PhiWpk=0.4d0     !Perry, et al 2015 RISR arc example
!!      PhiWpk=1d-5    !successful grad-drift attempts
!!!      PhiWpk=1d-4    !Swoboda blur testing
!!!      PhiWpk=0.05d0    !testing of convergent Hall drifts
!!!      PhiWpk=5d0
!!!      sigPhiWx3=100d3
!!!      meanPhiWx3=0d0
!!
!!!      W0pk=0.3d3
!!!      sigW0x3=100d3
!!!      meanW0x3=0d0
!!!      PhiWpk=2d0
!!!      sigPhiWx3=100d3
!!!      meanPhiWx3=0d0
!!
!!    sigx2=50d3
!!    meanx2=0d0
!!!    sigx3=10d3
!!    sigx3=25d3
!!    meant=900d0
!!    sigt=450d0
!!    x30amp=0d3
!!    varc=200d0
!!
!!    !DISTURBANCE ELECTRON PRECIPITATION PATTERN
!!    do ix3=1,lx3
!!      do ix2=1,lx2
!!        meanx3=-varc*meant+varc*t
!!
!!        W0(ix2,ix3,2)=W0pk
!!!        PhiWmWm2(ix2,ix3)=PhiWpk
!!        PhiWmWm2(ix2,ix3,2)=(PhiWpk*exp(-(x3(ix3)-(-sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x2(ix2)) ))**2/2d0/sigx3**2)- &
!!                          Phiwpk*exp(-(x3(ix3)-(sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x2(ix2)) ))**2/2d0/sigx3**2) )* &
!!                            exp(-(x2(ix2)-meanx2)**2/2d0/sigx2**2)*exp(-(t-meant)**2/2d0/sigt**2)
!!        if (PhiWmWm2(ix2,ix3,2) < 0d0) then
!!          PhiWmWm2(ix2,ix3,2)=0d0
!!        end if
!!!        PhiWmWm2(ix2,ix3,2)=(PhiWpk*exp(-(x3(ix3)-(-sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x2(ix2)) ))**2/2d0/sigx3**2) )* &
!!!                              exp(-(x2(ix2)-meanx2)**2/2d0/sigx2**2)*exp(-(t-meant)**2/2d0/sigt**2)
!!      end do
!!    end do
!
!  end subroutine precipBCs_original
!

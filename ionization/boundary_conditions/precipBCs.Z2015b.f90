  subroutine precipBCs(t,x2,x3,W0,PhiWmWm2)

    !------------------------------------------------------------
    !-------LOAD UP ARRAYS CONTAINING TOP BOUNDARY CHAR. ENERGY
    !-------AND TOTAL ENERGY FLUX.  GRID VARIABLES INCLUDE
    !-------GHOST CELLS
    !------------------------------------------------------------

    use phys_consts

    implicit none
    real(8), intent(in) :: t
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3
    real(8), dimension(:,:,:), intent(out) :: W0,PhiWmWm2

    real(8) :: W0pk,PhiWpk,meanW0x3,meanPhiWx3,sigW0x3,sigPhiWx3
    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve,sigt,meant
    integer :: ix2,ix3,iprec,lx2,lx3,lprec

    
    lx2=size(W0,1)
    lx3=size(W0,2)
    lprec=size(W0,3)    !assumed to be 2 in this subroutine


    !BACKGROUND PRECIPITATION
!    W0pk=5000d0
!    PhiWpk=0.01d0
    W0pk=3d3
    PhiWpk=1d-5
    do ix3=1,lx3
      do ix2=1,lx2
        W0(ix2,ix3,1)=W0pk 
        PhiWmWm2(ix2,ix3,1)=PhiWpk       
      end do
    end do


    !PARAMETERS FOR DISTURBANCE PRECIPITATION
    W0pk=100d0
    PhiWpk=0d0

    sigx2=50d3
    meanx2=0d0
    sigx3=25d3
    meant=900d0
    sigt=450d0
    x30amp=0d3
    varc=200d0


    !DISTURBANCE ELECTRON PRECIPITATION PATTERN
    do ix3=1,lx3
      do ix2=1,lx2
        meanx3=-varc*meant+varc*t

        W0(ix2,ix3,2)=W0pk
        PhiWmWm2(ix2,ix3,2)=(PhiWpk*exp(-(x3(ix3)-(-sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x2(ix2)) ))**2/2d0/sigx3**2)- &
                          Phiwpk*exp(-(x3(ix3)-(sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x2(ix2)) ))**2/2d0/sigx3**2) )* &
                            exp(-(x2(ix2)-meanx2)**2/2d0/sigx2**2)*exp(-(t-meant)**2/2d0/sigt**2)
        if (PhiWmWm2(ix2,ix3,2) < 0d0) then
          PhiWmWm2(ix2,ix3,2)=0d0
        end if
      end do
    end do

 end subroutine precipBCs

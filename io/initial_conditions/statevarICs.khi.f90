  subroutine statevarICs(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    implicit none
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x20amp,x30amp,varc,meanx2
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp
    real(8), dimension(1:size(x2,1)-4) :: x3dev
    real(8), dimension(1:size(x3all,1)-4) :: x2dev
    integer :: im
    integer, parameter :: lmodes=12    !number of cos modes to use to seed patch
    real(8) :: phase
    real(8), parameter :: pi=3.141592d0


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lsp=size(nsall,4)

    sigx2=8d3
    meanx3=10d3
    sigx3=10d3
    x30amp=0.1d3
    x20amp=0.1d3
    varc=0d0
    meanx2=0d0


    !SET UP NUMBER OF SEED MODES AND RANDOM PHASES
    call random_seed()   !default initialization (probably not really random but doesn't matter for now)
    x3dev=0d0
!    do im=1,lmodes    !comment out if no seed for instability
!      call random_number(phase)
!      phase=phase*2*pi
!      x3dev=x3dev+x30amp*cos(2d0*pi/(2d0*sigx2/real(im,8))*x2(1:lx2)+phase)   !center line for patch edge
!    end do

    call random_seed()   !default initialization (probably not really random but doesn't matter for now)
    x2dev=0d0
    do im=1,lmodes    !comment out if no seed for instability
      call random_number(phase)
      phase=phase*2*pi
      x2dev=x2dev+x20amp*cos(2d0*pi/(2d0*sigx3/real(im,8))*x3all(1:lx3all)+phase)   !center line for patch edge
    end do


    !CREATE THE PATCH - LEAVE TEMPERATURES AND PARALLEL FLOWS UNIFORM
    do isp=1,lsp
      do ix3=1,lx3all
        do ix2=1,lx2
          do ix1=1,lx1
            nsall(ix1,ix2,ix3,isp)=nsin(ix1,1,1,isp)+7.5d0*nsin(ix1,1,1,isp)* &
                                    exp(-1d0*(x3all(ix3)-meanx3-x3dev(ix2))**10/2d0/sigx3**10)* &
                                    exp(-1d0*(x2(ix2)-meanx2-x2dev(ix3))**10/2d0/sigx2**10)    !square patch
!            nsall(ix1,ix2,ix3,isp)=nsin(ix1,1,1,isp)+7.5d0*nsin(ix1,1,1,isp)* &
!                                    exp(-1d0*((x3all(ix3)-meanx3-x3dev(ix2))**2+(x2(ix2)-meanx2)**2)**5 &
!                                    /2d0/sigx3**10)    !round structure
          end do
        end do
      end do
    end do
    do isp=1,lsp
      vs1all(1:lx1,1:lx2,1:lx3all,isp)=spread(spread(vs1in(:,1,1,isp),2,lx2),3,lx3all)
    end do
    do isp=1,lsp
      Tsall(1:lx1,1:lx2,1:lx3all,isp)=spread(spread(Tsin(:,1,1,isp),2,lx2),3,lx3all)
    end do

  end subroutine statevarICs


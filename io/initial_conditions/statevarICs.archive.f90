  subroutine input_spread_copy(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lsp=size(nsall,4)

      do isp=1,lsp
        nsall(1:lx1,1:lx2,1:lx3all,isp)=spread(spread(nsin(:,1,1,isp),2,lx2),3,lx3all)
      end do
      do isp=1,lsp
        vs1all(1:lx1,1:lx2,1:lx3all,isp)=spread(spread(vs1in(:,1,1,isp),2,lx2),3,lx3all)
      end do
      do isp=1,lsp
        Tsall(1:lx1,1:lx2,1:lx3all,isp)=spread(spread(Tsin(:,1,1,isp),2,lx2),3,lx3all)
      end do

  end subroutine input_spread_copy


  subroutine input_spread_patch2(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp
    real(8), dimension(1:size(x2,1)-4) :: x3dev
    integer :: im
    integer, parameter :: lmodes=70    !number of cos modes to use to seed patch
    real(8) :: phase


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lsp=size(nsall,4)

    sigx2=12d3
!    meanx3=10d3
    meanx3=30d3
!    sigx3=10d3
    sigx3=20d3
    x30amp=0.5d3
!    x30amp=0.1d3
    varc=0d0
    meanx2=0d0
    x2enve=16d3


    !SET UP NUMBER OF SEED MODES AND RANDOM PHASES
    call random_seed()   !default initialization (probably not really random but doesn't matter for now)
    x3dev=0d0


    !CREATE THE PATCH - LEAVE TEMPERATURES AND PARALLEL FLOWS UNIFORM
    do im=1,lmodes
      call random_number(phase)   !single random phase for each mode
      phase=phase*2*pi

      do isp=1,lsp
        do ix3=1,lx3all
          do ix2=1,lx2
            do ix1=1,lx1
              !PUT IN THE PATCH
              nsall(ix1,ix2,ix3,isp)=nsin(ix1,1,1,isp)+nsin(ix1,1,1,isp)* &
                                    exp(-1d0*(x3all(ix3)-meanx3-x3dev(ix2))**10/2d0/sigx3**10)* &
                                    exp(-(x2(ix2)-meanx2)**6/2d0/x2enve**6)

              !ADD SOME PERTURBATION (10 pct. of local density)
              nsall(ix1,ix2,ix3,isp)=nsall(ix1,ix2,ix3,isp)+0.1* &
                (2d0*pi/2d0/sigx2)*(2d0*sigx2/real(im,8)/2d0/pi)*nsall(ix1,ix2,ix3,isp)* &
                sin(2d0*pi/(2d0*sigx2/real(im,8))*x2(ix2)+phase)
            end do
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

  end subroutine input_spread_patch2


  subroutine input_spread_patch3(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp
    real(8), dimension(1:size(x2,1)-4) :: x3dev
    integer :: im
    integer, parameter :: lmodes=12    !number of cos modes to use to seed patch
    real(8) :: phase


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lsp=size(nsall,4)

!    sigx2=12d3
    sigx2=10d3
    meanx3=10d3
!    meanx3=40d3
    sigx3=10d3
!    x30amp=1d3
    x30amp=0.1d3
    varc=0d0
    meanx2=0d0
    x2enve=16d3   !why not just sigx2?
!    x2enve=8d3


    !SET UP NUMBER OF SEED MODES AND RANDOM PHASES
    call random_seed()   !default initialization (probably not really random but doesn't matter for now)
    x3dev=0d0
    do im=1,lmodes    !comment out if no seed for instability
      call random_number(phase)
      phase=phase*2*pi
      x3dev=x3dev+x30amp*cos(2d0*pi/(2d0*sigx2/real(im,8))*x2(1:lx2)+phase)   !center line for patch edge
    end do


    !CREATE THE PATCH - LEAVE TEMPERATURES AND PARALLEL FLOWS UNIFORM
    do isp=1,lsp
      do ix3=1,lx3all
        do ix2=1,lx2
          do ix1=1,lx1
            nsall(ix1,ix2,ix3,isp)=nsin(ix1,1,1,isp)+7.5d0*nsin(ix1,1,1,isp)* &
                                    exp(-1d0*(x3all(ix3)-meanx3-x3dev(ix2))**10/2d0/sigx3**10)* &
                                    exp(-1d0*(x2(ix2)-meanx2)**6/2d0/x2enve**6)
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

  end subroutine input_spread_patch3


  subroutine input_spread_patch(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp
    real(8), dimension(1:size(x2,1)-4) :: x3dev
    integer :: im
    integer, parameter :: lmodes=12    !number of cos modes to use to seed patch
    real(8) :: phase


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lsp=size(nsall,4)

!    sigx2=12d3
    sigx2=10d3
    meanx3=10d3
!    meanx3=40d3
    sigx3=10d3
!    x30amp=1d3
    x30amp=0.1d3
    varc=0d0
    meanx2=0d0
    x2enve=16d3   !why not just sigx2?


    !SET UP NUMBER OF SEED MODES AND RANDOM PHASES
    call random_seed()   !default initialization (probably not really random but doesn't matter for now)
    x3dev=0d0
    do im=1,lmodes    !comment out if no seed for instability
      call random_number(phase)
      phase=phase*2*pi
      x3dev=x3dev+x30amp*cos(2d0*pi/(2d0*sigx2/real(im,8))*x2(1:lx2)+phase)   !center line for patch edge
    end do


    !CREATE THE PATCH - LEAVE TEMPERATURES AND PARALLEL FLOWS UNIFORM
    do isp=1,lsp
      do ix3=1,lx3all
        do ix2=1,lx2
          do ix1=1,lx1
            nsall(ix1,ix2,ix3,isp)=nsin(ix1,1,1,isp)+7.5d0*nsin(ix1,1,1,isp)* &
                                    exp(-1d0*(x3all(ix3)-meanx3-x3dev(ix2))**10/2d0/sigx3**10)* &
                                    exp(-1d0*(x2(ix2)-meanx2)**6/2d0/x2enve**6)
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

  end subroutine input_spread_patch


  subroutine input_spread_patch_khi(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

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

  end subroutine input_spread_patch_khi


  subroutine input_spread_edge(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lsp=size(nsall,4)

!    sigx2=60d3
!    meanx3=175d3
!!    sigx3=20d3
!    sigx3=2.5d3
!    x30amp=10d3
!    varc=0d0
!    meanx2=0d0
!    x2enve=80d3

    sigx2=12d3
    meanx3=17.5d3
    sigx3=250d0
    x30amp=1d3
    varc=0d0
    meanx2=0d0
    x2enve=16d3


    do isp=1,lsp
      do ix3=1,lx3all
        do ix2=1,lx2
          do ix1=1,lx1
            if (x1(ix1)>190d3) then
!              nsall(ix1,ix2,ix3,isp)=nsin(ix1,1,1,isp)+10*nsin(ix1,1,1,isp)*1d0/2d0-10*nsin(ix1,1,1,isp)*1d0/2d0* &
!                                      erf((x3all(ix3)-meanx3-x30amp*cos(2d0*pi/(sigx2)*x2(ix2))) &
!                                      /sqrt(2d0)/sigx3)*exp(-(x2(ix2)-meanx2)**6/2.0/x2enve**6)
              nsall(ix1,ix2,ix3,isp)=nsin(ix1,1,1,isp)+3d0*nsin(ix1,1,1,isp)*1d0/2d0-3d0*nsin(ix1,1,1,isp)*1d0/2d0* &
                                      erf((x3all(ix3)-meanx3-x30amp*cos(2d0*pi/(sigx2)*x2(ix2))) &
                                      /sqrt(2d0)/sigx3)*exp(-(x2(ix2)-meanx2)**6/2.0/x2enve**6)
            else    !chop out lower alts. to prevent E-region shorting...
              nsall(ix1,ix2,ix3,isp)=1d4
            end if
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


  end subroutine input_spread_edge


  subroutine input_spread_extend(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp
    integer :: lx3in


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lx3in=size(nsin,3)
    lsp=size(nsall,4)


    !LX1,LX2 ASSUMED TO BE SAME
    if (lx3all >= lx3in) then
      do isp=1,lsp
        nsall(1:lx1,1:lx2,1:lx3in,isp)=nsin(:,:,1:lx3in,isp)
      end do
      do isp=1,lsp
        vs1all(1:lx1,1:lx2,1:lx3in,isp)=vs1in(:,:,1:lx3in,isp)
      end do
      do isp=1,lsp
        Tsall(1:lx1,1:lx2,1:lx3in,isp)=Tsin(:,:,1:lx3in,isp)
      end do

      do isp=1,lsp
        do ix3=lx3in+1,lx3all    !just repeat last entry for extras
          nsall(1:lx1,1:lx2,ix3,isp)=nsin(:,:,lx3in,isp)
          vs1all(1:lx1,1:lx2,ix3,isp)=vs1in(:,:,lx3in,isp)
          Tsall(1:lx1,1:lx2,ix3,isp)=Tsin(:,:,lx3in,isp)
        end do
      end do
    else
      do isp=1,lsp
        nsall(1:lx1,1:lx2,1:lx3all,isp)=nsin(:,:,1:lx3all,isp)
      end do
      do isp=1,lsp
        vs1all(1:lx1,1:lx2,1:lx3all,isp)=vs1in(:,:,1:lx3all,isp)
      end do
      do isp=1,lsp
        Tsall(1:lx1,1:lx2,1:lx3all,isp)=Tsin(:,:,1:lx3all,isp)
      end do
    end if

  end subroutine input_spread_extend


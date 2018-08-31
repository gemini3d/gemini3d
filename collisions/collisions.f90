module collisions

use phys_consts
implicit none

!THESE DATA BLOCKS ARE UGLY; COULD THERE BE A BETTER WAY?
real(8), dimension(lsp,ln) :: Csn
data Csn(1,:) /-1.0, 6.82d-10, 6.64d-10, -1.0/
data Csn(2,:) /2.44d-10, 4.34d-10, 4.27d-10, 0.69d-10/
data Csn(3,:) /2.58d-10, -1.0, 4.49d-10, 0.74d-10/
data Csn(4,:) /2.31d-10, 4.13d-10, -1.0, 0.65d-10/
data Csn(5,:) /4.42d-10, 7.47d-10, 7.25d-10, 1.45d-10/
data Csn(6,:) /-1.0, 33.6d-10, 32.0d-10, -1.0/
data Csn(7,:) /-1.0, -1.0, -1.0, -1.0/

real(8), dimension(lsp,ln) :: C2sn1
data C2sn1(1,:) /3.67d-11, 0.0, 0.0, 4.63d-12/
data C2sn1(2,:) /0.0, 0.0, 0.0, 0.0/
data C2sn1(3,:) /0.0, 5.14d-11, 0.0, 0.0/
data C2sn1(4,:) /0.0, 0.0, 2.59d-11, 0.0/
data C2sn1(5,:) /0.0, 0.0, 0.0, 0.0/
data C2sn1(6,:) /6.61d-11, 0.0, 0.0, 2.65d-10/
data C2sn1(7,:) /-1.0, -1.0, -1.0, -1.0/

real(8), dimension(lsp,ln) :: C2sn2
data C2sn2(1,:) /0.064, 0.0, 0.0, -1.0/
data C2sn2(2,:) /0.0, 0.0, 0.0, 0.0/
data C2sn2(3,:) /0.0, 0.069, 0.0, 0.0/
data C2sn2(4,:) /0.0, 0.0, 0.073, 0.0/
data C2sn2(5,:) /0.0, 0.0, 0.0, 0.0/
data C2sn2(6,:) /0.047, 0.0, 0.0, 0.083/
data C2sn2(7,:) /-1.0, -1.0, -1.0, -1.0/

real(8), dimension(lsp,lsp) :: Csj
data Csj(1,:) /0.22, 0.26, 0.25, 0.26, 0.22, 0.077, 1.87e-3/
data Csj(2,:) /0.14, 0.16, 0.16, 0.17, 0.13, 0.042, 9.97e-4/
data Csj(3,:) /0.15, 0.17, 0.17, 0.18, 0.14, 0.045, 1.07e-3/
data Csj(4,:) /0.13, 0.16, 0.15, 0.16, 0.12, 0.039, 9.347e-4/
data Csj(5,:) /0.25, 0.28, 0.28, 0.28, 0.24, 0.088, 2.136e-3/
data Csj(6,:) /1.23, 1.25, 1.25, 1.25, 1.23, 0.90,  29.7e-3/
data Csj(7,:) /54.5, 54.5, 54.5, 54.5, 54.5, 54.5,  38.537/

contains


  subroutine maxwell_colln(isp,isp2,nn,Tn,Ts,nusn)

    !------------------------------------------------------------
    !-------COMPUTE MAXWELL COLLISIONS OF ISP WITH ISP2.  ION 
    !-------TEMPERATURE/DENSITY ARRAYS EXPECTED TO INCLUDE GHOST CELLS
    !------------------------------------------------------------
    !-------Note that it is done on a per species basis

    integer, intent(in) :: isp,isp2
    real(8), dimension(:,:,:,:), intent(in) :: nn
    real(8), dimension(:,:,:), intent(in) :: Tn
    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: Ts

    real(8), dimension(1:size(Tn,1),1:size(Tn,2),1:size(Tn,3)), &
              intent(out) :: nusn

    integer :: lx1,lx2,lx3
    real(8) :: mred
    real(8),dimension(1:size(Tn,1),1:size(Tn,2),1:size(Tn,3)) :: Teff

    lx1=size(Ts,1)-4
    lx2=size(Ts,2)-4
    lx3=size(Ts,3)-4


    if (isp<lsp) then    !ion-neutral
      if (Csn(isp,isp2)<0.0) then    !resonant
        if (isp==1 .and. isp2==4) then
          Teff=Tn+Ts(1:lx1,1:lx2,1:lx3,isp)/16.0
          nusn=C2sn1(isp,isp2)*Teff**0.5*nn(:,:,:,isp2)*1d-6
        else
          Teff=0.5*(Tn+Ts(1:lx1,1:lx2,1:lx3,isp))
          nusn=C2sn1(isp,isp2)*(1.0-C2sn2(isp,isp2)*log10(Teff))**2.0* &
               (Teff**0.5)*nn(:,:,:,isp2)*1d-6
        end if
      else    !nonresonant
        nusn=Csn(isp,isp2)*nn(:,:,:,isp2)*1d-6
      end if
    else    !electron-neutral
      Teff=Ts(1:lx1,1:lx2,1:lx3,isp)

      select case (isp2)
        case (1)
          nusn=8.9d-11*(1.0+5.7d-4*Teff)*(Teff**0.5)*nn(:,:,:,isp2)*1d-6
        case (2)
          nusn=2.33d-11*(1.0-1.21d-4*Teff)*(Teff)*nn(:,:,:,isp2)*1d-6
        case (3)
          nusn=1.82d-10*(1.0+3.6d-2*(Teff**0.5))*(Teff**0.5)*nn(:,:,:,isp2)*1d-6
        case (4)
          nusn=4.5d-9*(1.0-1.35d-4*Teff)*(Teff**0.5)*nn(:,:,:,isp2)*1d-6
      end select
    end if

  end subroutine maxwell_colln


  subroutine coulomb_colln(isp,isp2,ns,Ts,vs1,nusj,Phisj,Psisj)

    !------------------------------------------------------------
    !-------COMPUTE COULOMB COLLISIONS OF ISP WITH ISP2.  
    !-------TEMPERATURE/DENSITY ARRAYS EXPECTED TO INCLUDE GHOST CELLS
    !-------NOTE THAT OTHER PIECES OF THE CODE REQUIRE SELF COLLISIONS
    !-------TO BE ZERO TO YIELD CORRECT OUTPUT (SOURCES.MOD)
    !------------------------------------------------------------
    !-------Note that it is done on a per species basis

    integer, intent(in) :: isp,isp2
    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1

    real(8), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), &
              intent(out) :: nusj,Phisj,Psisj

    integer :: lx1,lx2,lx3
    real(8) :: mred
    real(8),dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4) &
              :: Teff,Wsj,Phitmp


    lx1=size(Ts,1)-4
    lx2=size(Ts,2)-4
    lx3=size(Ts,3)-4

    if (isp==isp2) then     !zero out all self collision terms (would need to be changed if non-Maxwellian distribution used).
      nusj=0d0
      Phisj=0d0
      Psisj=0d0
    else
      Teff=(ms(isp2)*Ts(1:lx1,1:lx2,1:lx3,isp)+ms(isp)* &
             Ts(1:lx1,1:lx2,1:lx3,isp2))/(ms(isp2)+ms(isp))
      nusj=Csj(isp,isp2)*ns(1:lx1,1:lx2,1:lx3,isp2)*1d-6/Teff**1.5d0

      mred=ms(isp)*ms(isp2)/(ms(isp)+ms(isp2))
      Wsj=abs(vs1(1:lx1,1:lx2,1:lx3,isp)-vs1(1:lx1,1:lx2,1:lx3,isp2))/ &
            sqrt(2*kB*Teff/mred)
      Psisj=exp(-Wsj**2d0)
      where (Wsj<0.1d0)
        Phisj=1d0
      elsewhere
        Phisj=3d0/4d0*sqrt(pi)*erf(Wsj)/Wsj**3d0-3d0/2d0/Wsj**2d0*Psisj
      end where
    end if
  end subroutine coulomb_colln


  subroutine thermal_conduct(isp,Ts,ns,nn,J1,lambda,beta)

    !------------------------------------------------------------
    !-------COMPUTE THERMAL CONDUCTIVITY.  TEMPERATURE ARRAY
    !-------IS EXPECTED TO INCLUDE GHOST CELLS
    !------------------------------------------------------------
    !-------Note that it is done on a per species basis

    integer, intent(in) :: isp
    real(8), dimension(-1:,-1:,-1:), intent(in) :: Ts,ns

    real(8), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4,ln), intent(in) :: nn
    real(8), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), intent(in) :: J1

    real(8), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), intent(out) :: lambda,beta

    integer :: lx1,lx2,lx3


    lx1=size(Ts,1)-4
    lx2=size(Ts,2)-4
    lx3=size(Ts,3)-4

    if (isp<lsp) then       !ion species
      lambda=25.0/8.0*kB**2*Ts(1:lx1,1:lx2,1:lx3)**(5.0/2.0)/ms(isp)/(Csj(isp,isp)*1d-6)
      beta=0.0
    else                  !electrons
      lambda=elchrg*100.0*7.7d5*Ts(1:lx1,1:lx2,1:lx3)**(5.0/2.0)/ &
          (1.0+3.22d4*Ts(1:lx1,1:lx2,1:lx3)**2/ns(1:lx1,1:lx2,1:lx3)* &
            (nn(:,:,:,1)*1.1d-16*(1+5.7d-4*Ts(1:lx1,1:lx2,1:lx3)) + &
            nn(:,:,:,2)*2.82d-17*sqrt(Ts(1:lx1,1:lx2,1:lx3))* &
            (1-1.21d-4*Ts(1:lx1,1:lx2,1:lx3))+nn(:,:,:,3)* &
            2.2d-16*(1+3.6d-2*sqrt(Ts(1:lx1,1:lx2,1:lx3))) ))
      beta=5.0/2.0*kB/elchrg*J1  
    end if

  end subroutine thermal_conduct


  subroutine conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,muPvn,muHvn) 

    !------------------------------------------------------------
    !-------COMPUTE THE CONDUCTIVITIES OF THE IONOSPHERE.  STATE
    !-------VARS. INCLUDE GHOST CELLS
    !------------------------------------------------------------

    real(8), dimension(:,:,:,:), intent(in) :: nn
    real(8), dimension(:,:,:), intent(in) :: Tn
    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
    real(8), dimension(-1:,-1:,-1:), intent(in) :: B1
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), intent(out) :: sig0,sigP,sigH
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,lsp), intent(out) :: muP,muH,muPvn,muHvn

    integer :: isp,isp2,lx1,lx2,lx3
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: OMs
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: nu,nuej,Phisj,Psisj,nutmp,mupar,mubase,rho


    lx1=size(Ts,1)-4
    lx2=size(Ts,2)-4
    lx3=size(Ts,3)-4


    !MOBILITIES
    do isp=1,lsp
!      OMs=qs(isp)*abs(B1)/ms(isp)                 !cyclotron, abs() is sketch, needs to be checked.  Basically a negative sign here is fine, while abs messes up direction of Hall current
      OMs=qs(isp)*B1(1:lx1,1:lx2,1:lx3)/ms(isp)                 !cyclotron, a negative sign from B1 here is fine for cartesian, but for dipole this should be the magnitude since the magnetic field is *assumed* to be along the x1-direction
   
      nu=0d0
      do isp2=1,ln
        call maxwell_colln(isp,isp2,nn,Tn,Ts,nutmp)
        nu=nu+nutmp
      end do

      if (isp<lsp) then
        mubase=qs(isp)/ms(isp)/nu      !parallel mobility
      else
        nuej=0d0
        do isp2=1,lsp
          call coulomb_colln(isp,isp2,ns,Ts,vs1,nutmp,Phisj,Psisj)
          nuej=nuej+nutmp
        end do

        mupar=qs(lsp)/ms(lsp)/(nu+nuej)
        mubase=qs(lsp)/ms(lsp)/nu
      end if

      !modified mobilities for neutral wind calculations
      muPvn(:,:,:,isp)=nu**2/(nu**2+OMs**2)
      muHvn(:,:,:,isp)=-1d0*nu*OMs/(nu**2+OMs**2)

      !full mobilities   
      muP(:,:,:,isp)=mubase*muPvn(:,:,:,isp)           !Pederson
      muH(:,:,:,isp)=mubase*muHvn(:,:,:,isp)       !Hall
    end do


    !CONDUCTIVITIES
    sig0=ns(1:lx1,1:lx2,1:lx3,lsp)*qs(lsp)*mupar    !parallel includes only electrons...

    sigP=0d0
    sigH=0d0
    do isp=1,lsp
      rho=ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)
      sigP=sigP+rho*muP(:,:,:,isp)
      sigH=sigH+rho*muH(:,:,:,isp)
    end do

!    sigH=max(sigH,0d0)    !to deal with precision issues.  This actually causes errors in Cartesian northern hemisphere grids...
  end subroutine conductivities


  subroutine capacitance(ns,B1,flagcap,incap) 

    !------------------------------------------------------------
    !-------COMPUTE THE INERTIAL CAPACITANCE OF THE IONOSPHERE.
    !-------DENSITY/MAG FIELD STATE VARIABLE INCLUDES GHOST CELLS.
    !------------------------------------------------------------

    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: ns
    real(8), dimension(-1:,-1:,-1:), intent(in) :: B1
    integer, intent(in) :: flagcap

    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), intent(out) :: incap

    integer :: lx1,lx2,lx3,isp


    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4

    incap=0d0
    do isp=1,lsp
      incap=incap+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)
    end do

    incap=incap/B1(1:lx1,1:lx2,1:lx3)**2

    if (flagcap==2) then
      write(*,*) '!!! Augmenting capacitance with a magnetospheric contribution...'
      incap=incap+30d0/980d3    !kludge the value to account for a magnetosheric contribution - this is just a random guess that makes the KHI examples work well; a better value should be investigated
    end if
  end subroutine capacitance


end module collisions

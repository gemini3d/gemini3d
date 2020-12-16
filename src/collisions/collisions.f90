module collisions

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp, lsp, ln, ms, kb, pi, elchrg, qs, debug
use config, only: gemini_cfg

implicit none (type, external)
private
public :: thermal_conduct, conductivities, capacitance, maxwell_colln, coulomb_colln


real(wp), parameter :: Csn(lsp,ln) = reshape( &
[-1.0_wp, 6.82e-10_wp, 6.64e-10_wp, -1.0_wp, &
2.44e-10_wp, 4.34e-10_wp, 4.27e-10_wp, 0.69e-10_wp, &
2.58e-10_wp, -1.0_wp, 4.49e-10_wp, 0.74e-10_wp, &
2.31e-10_wp, 4.13e-10_wp, -1.0_wp, 0.65e-10_wp, &
4.42e-10_wp, 7.47e-10_wp, 7.25e-10_wp, 1.45e-10_wp, &
-1.0_wp, 33.6e-10_wp, 32.0e-10_wp, -1.0_wp, &
-1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp], shape(Csn), order=[2,1])

real(wp), parameter :: C2sn1(lsp,ln) = reshape( &
[3.67e-11_wp, 0._wp, 0._wp, 4.63e-12_wp, &
0._wp, 0._wp, 0._wp, 0._wp, &
0._wp, 5.14e-11_wp, 0._wp, 0._wp, &
0._wp, 0._wp, 2.59e-11_wp, 0._wp, &
0._wp, 0._wp, 0._wp, 0._wp, &
6.61e-11_wp, 0._wp, 0._wp, 2.65e-10_wp, &
-1._wp, -1._wp, -1._wp, -1._wp], shape(C2sn1), order=[2,1])

real(wp), parameter :: C2sn2(lsp,ln) = reshape( &
[0.064_wp, 0._wp, 0._wp, -1._wp, &
0._wp, 0._wp, 0._wp, 0._wp, &
0._wp, 0.069_wp, 0._wp, 0._wp, &
0._wp, 0._wp, 0.073_wp, 0._wp, &
0._wp, 0._wp, 0._wp, 0._wp, &
0.047_wp, 0._wp, 0._wp, 0.083_wp, &
-1._wp, -1._wp, -1._wp, -1._wp], shape(C2sn2), order=[2,1])

real(wp), parameter :: Csj(lsp,lsp) = reshape( &
[0.22_wp, 0.26_wp, 0.25_wp, 0.26_wp, 0.22_wp, 0.077_wp, 1.87e-3_wp, &
0.14_wp, 0.16_wp, 0.16_wp, 0.17_wp, 0.13_wp, 0.042_wp, 9.97e-4_wp, &
0.15_wp, 0.17_wp, 0.17_wp, 0.18_wp, 0.14_wp, 0.045_wp, 1.07e-3_wp, &
0.13_wp, 0.16_wp, 0.15_wp, 0.16_wp, 0.12_wp, 0.039_wp, 9.347e-4_wp, &
0.25_wp, 0.28_wp, 0.28_wp, 0.28_wp, 0.24_wp, 0.088_wp, 2.136e-3_wp, &
1.23_wp, 1.25_wp, 1.25_wp, 1.25_wp, 1.23_wp, 0.90_wp,  29.7e-3_wp, &
54.5_wp, 54.5_wp, 54.5_wp, 54.5_wp, 54.5_wp, 54.5_wp,  38.537_wp], shape(Csj), order=[2,1])

real(wp), parameter :: thermal_coeff(lsp-1) = [0.1019e-12_wp,0.0747e-12_wp,0.0754e-12_wp,0.0701e-12_wp, &
                                        0.1068e-12_wp,0.3986e-12_wp]

contains


subroutine maxwell_colln(isp,isp2,nn,Tn,Ts,nusn)

!------------------------------------------------------------
!-------COMPUTE MAXWELL COLLISIONS OF ISP WITH ISP2.  ION
!-------TEMPERATURE/DENSITY ARRAYS EXPECTED TO INCLUDE GHOST CELLS
!------------------------------------------------------------
!-------Note that it is done on a per species basis

integer, intent(in) :: isp,isp2
real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: Ts

real(wp), dimension(1:size(Tn,1),1:size(Tn,2),1:size(Tn,3)), &
          intent(out) :: nusn

integer :: lx1,lx2,lx3
real(wp) :: mred
real(wp),dimension(1:size(Tn,1),1:size(Tn,2),1:size(Tn,3)) :: Teff

lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4


if (isp<lsp) then    !ion-neutral
  if (Csn(isp,isp2)<0.0) then    !resonant
    if (isp==1 .and. isp2==4) then
      Teff=Tn+Ts(1:lx1,1:lx2,1:lx3,isp)/16.0
      nusn=C2sn1(isp,isp2)*Teff**0.5*nn(:,:,:,isp2)*1e-6_wp
    else
      Teff=0.5*(Tn+Ts(1:lx1,1:lx2,1:lx3,isp))
      nusn=C2sn1(isp,isp2)*(1.0-C2sn2(isp,isp2)*log10(Teff))**2.0* &
           (Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
    end if
  else    !nonresonant
    nusn=Csn(isp,isp2)*nn(:,:,:,isp2)*1e-6_wp
  end if
else    !electron-neutral
  Teff=Ts(1:lx1,1:lx2,1:lx3,isp)

  select case (isp2)
    case (1)
      nusn=8.9e-11_wp*(1.0+5.7e-4_wp*Teff)*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
    case (2)
      nusn=2.33e-11_wp*(1.0-1.21e-4_wp*Teff)*(Teff)*nn(:,:,:,isp2)*1e-6_wp
    case (3)
      nusn=1.82e-10_wp*(1.0+3.6e-2_wp*(Teff**0.5))*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
    case (4)
      nusn=4.5e-9_wp*(1.0-1.35e-4_wp*Teff)*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
    case default
      write(stderr,*) 'ERROR: isp2 value is unknown: ',isp2
      error stop
  end select
end if

end subroutine maxwell_colln


pure subroutine coulomb_colln(isp,isp2,ns,Ts,vs1,nusj,Phisj,Psisj)

!------------------------------------------------------------
!-------COMPUTE COULOMB COLLISIONS OF ISP WITH ISP2.
!-------TEMPERATURE/DENSITY ARRAYS EXPECTED TO INCLUDE GHOST CELLS
!-------NOTE THAT OTHER PIECES OF THE CODE REQUIRE SELF COLLISIONS
!-------TO BE ZERO TO YIELD CORRECT OUTPUT (SOURCES.MOD)
!------------------------------------------------------------
!-------Note that it is done on a per species basis

integer, intent(in) :: isp,isp2
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1

real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), &
          intent(out) :: nusj,Phisj,Psisj

integer :: lx1,lx2,lx3
real(wp) :: mred
real(wp),dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4) &
          :: Teff,Wsj,Phitmp


lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4

if (isp==isp2) then     !zero out all self collision terms (would need to be changed if non-Maxwellian distribution used).
  nusj = 0
  Phisj = 0
  Psisj = 0
else
  Teff=(ms(isp2)*Ts(1:lx1,1:lx2,1:lx3,isp)+ms(isp)* &
         Ts(1:lx1,1:lx2,1:lx3,isp2))/(ms(isp2)+ms(isp))
  nusj=Csj(isp,isp2)*ns(1:lx1,1:lx2,1:lx3,isp2)*1e-6_wp/Teff**1.5_wp

  mred=ms(isp)*ms(isp2)/(ms(isp)+ms(isp2))
  Wsj=abs(vs1(1:lx1,1:lx2,1:lx3,isp)-vs1(1:lx1,1:lx2,1:lx3,isp2))/ &
        sqrt(2*kB*Teff/mred)
  Psisj=exp(-Wsj**2.0_wp)
  where (Wsj<0.1_wp)
    Phisj=1.0_wp
  elsewhere
    Phisj=3.0_wp/4.0_wp*sqrt(pi)*erf(Wsj)/Wsj**3.0_wp-3.0_wp/2.0_wp/Wsj**2.0_wp*Psisj
  end where
end if
end subroutine coulomb_colln


pure subroutine thermal_conduct(isp,Ts,ns,nn,J1,lambda,beta)

!------------------------------------------------------------
!-------COMPUTE THERMAL CONDUCTIVITY.  TEMPERATURE ARRAY
!-------IS EXPECTED TO INCLUDE GHOST CELLS
!------------------------------------------------------------
!-------Note that it is done on a per species basis

integer, intent(in) :: isp
real(wp), dimension(-1:,-1:,-1:), intent(in) :: Ts,ns

real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4,ln), intent(in) :: nn
real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), intent(in) :: J1

real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), intent(out) :: lambda,beta

integer :: lx1,lx2,lx3


lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4

if (isp<lsp) then       !ion species
!  lambda=25.0_wp/8.0_wp*kB**2*Ts(1:lx1,1:lx2,1:lx3)**(5.0_wp/2.0_wp)/ms(isp)/(Csj(isp,isp)*1e-6_wp)
  ! avoids precision issues by precomputing the transport coefficients (see parameter blocks above)
  lambda=thermal_coeff(isp)*Ts(1:lx1,1:lx2,1:lx3)**(5.0_wp/2.0_wp)
  beta=0.0
else                  !electrons
  lambda=elchrg*100.0_wp*7.7e5_wp*Ts(1:lx1,1:lx2,1:lx3)**(5.0_wp/2.0_wp)/ &
      (1.0+3.22e4_wp*Ts(1:lx1,1:lx2,1:lx3)**2/ns(1:lx1,1:lx2,1:lx3)* &
        (nn(:,:,:,1)*1.1e-16_wp*(1+5.7e-4_wp*Ts(1:lx1,1:lx2,1:lx3)) + &
        nn(:,:,:,2)*2.82e-17_wp*sqrt(Ts(1:lx1,1:lx2,1:lx3))* &
        (1-1.21e-4_wp*Ts(1:lx1,1:lx2,1:lx3))+nn(:,:,:,3)* &
        2.2e-16_wp*(1+3.6e-2_wp*sqrt(Ts(1:lx1,1:lx2,1:lx3))) ))
  beta=5.0_wp/2.0_wp*kB/elchrg*J1
end if

end subroutine thermal_conduct


subroutine conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,muPvn,muHvn,sigPgrav,sigHgrav)

!------------------------------------------------------------
!-------COMPUTE THE CONDUCTIVITIES OF THE IONOSPHERE.  STATE
!-------VARS. INCLUDE GHOST CELLS
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), intent(out) :: sig0,sigP,sigH
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,lsp), intent(out) :: muP,muH,muPvn,muHvn

integer :: isp,isp2,lx1,lx2,lx3
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: OMs
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: nu,nuej,Phisj,Psisj,nutmp,mupar,mubase,rho
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), intent(out) :: sigPgrav,sigHgrav

lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4


!MOBILITIES
do isp=1,lsp
!      OMs=qs(isp)*abs(B1)/ms(isp)
  !! cyclotron, abs() is sketch, needs to be checked.
  !! Basically a negative sign here is fine, while abs messes up direction of Hall current
  OMs=qs(isp)*B1(1:lx1,1:lx2,1:lx3)/ms(isp)
  !! cyclotron, a negative sign from B1 here is fine for cartesian, but for dipole this should be the magnitude
  !! since the magnetic field is *assumed* to be along the x1-direction

  nu = 0
  do isp2=1,ln
    call maxwell_colln(isp,isp2,nn,Tn,Ts,nutmp)
    nu=nu+nutmp
  end do

  if (isp<lsp) then
    mubase=qs(isp)/ms(isp)/nu      !parallel mobility
  else
    nuej = 0
    do isp2=1,lsp
      call coulomb_colln(isp,isp2,ns,Ts,vs1,nutmp,Phisj,Psisj)
      nuej=nuej+nutmp
    end do

    mupar=qs(lsp)/ms(lsp)/(nu+nuej)
    mubase=qs(lsp)/ms(lsp)/nu
  end if

  !modified mobilities for neutral wind calculations
  muPvn(:,:,:,isp)=nu**2/(nu**2+OMs**2)
  muHvn(:,:,:,isp)=-1.0_wp*nu*OMs/(nu**2+OMs**2)

  !electrical mobilities
  muP(:,:,:,isp)=mubase*muPvn(:,:,:,isp)           !Pederson
  muH(:,:,:,isp)=mubase*muHvn(:,:,:,isp)       !Hall

  !gravity mobilities???
end do


!CONDUCTIVITIES
sig0=ns(1:lx1,1:lx2,1:lx3,lsp)*qs(lsp)*mupar    !parallel includes only electrons...

sigP = 0
sigH = 0
do isp=1,lsp
  rho=ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)    !rho is charge density here
  sigP=sigP+rho*muP(:,:,:,isp)
  sigH=sigH+rho*muH(:,:,:,isp)
end do

!    sigH=max(sigH,0.0_wp)    !to deal with precision issues.  This actually causes errors in Cartesian northern hemisphere grids...


!Gravitational "conductivities"
sigPgrav=0.0_wp
sigHgrav=0.0_wp
do isp=1,lsp
  rho=ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)     !here, rho is used as mass density
  sigPgrav=sigPgrav+rho*muP(:,:,:,isp)
  sigHgrav=sigHgrav+rho*muH(:,:,:,isp)
end do

end subroutine conductivities


subroutine capacitance(ns,B1,cfg,incap)

!------------------------------------------------------------
!-------COMPUTE THE INERTIAL CAPACITANCE OF THE IONOSPHERE.
!-------DENSITY/MAG FIELD STATE VARIABLE INCLUDES GHOST CELLS.
!------------------------------------------------------------

real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns
real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
type(gemini_cfg), intent(in) :: cfg
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), intent(out) :: incap

integer :: lx1,lx2,lx3,isp


lx1=size(ns,1)-4
lx2=size(ns,2)-4
lx3=size(ns,3)-4

incap = 0
do isp=1,lsp
  incap=incap+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)
end do

incap=incap/B1(1:lx1,1:lx2,1:lx3)**2

if (cfg%flagcap==2) then
  if (debug) print *, '!!! Augmenting capacitance with a magnetospheric contribution...'
  incap=incap + cfg%magcap / 980e3_wp
  !! augment the value to account for a magnetosheric contribution, based on user input.  Probably should
  !! be in the 5-35 F range...  Note that the assumes that the grid extends from ~90-1000 km in altitude approx.  
end if
end subroutine capacitance


end module collisions

module sources

use calculus, only : grad3d1
use collisions, only:  maxwell_colln, coulomb_colln
use phys_consts, only: wp, lsp, amu, kb, qs, ln, ms, gammas, elchrg, mn
use meshobj, only : curvmesh

implicit none (type, external)
private
public :: srcsenergy, srcsmomentum, srcscontinuity

interface srcsMomentum
  module procedure srcsMomentum_curv
end interface srcsMomentum

contains


pure subroutine srcsContinuity(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts, Pr, Lo)

!------------------------------------------------------------
!-------POPULATE SOURCE/LOSS ARRAYS FOR CONTINUITY EQUATION.  ION
!-------PARAMETER ARGUMENTS (AND GRID STUFF) SHOULD INCLUDE GHOST CELLS
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(inout) :: Pr,Lo
!! intent(out)

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: betanow,kreac,Teff,Te,dv2
integer :: lx1,lx2,lx3

lx1=size(ns,1)-4
lx2=size(ns,2)-4
lx3=size(ns,3)-4

Pr=0
Lo=0
Te=Ts(1:lx1,1:lx2,1:lx3,lsp)    !< Used in calculation of Lo
dv2=(vs1(1:lx1,1:lx2,1:lx3,1)-vn1)**2+(vs2(1:lx1,1:lx2,1:lx3,1)-vn2)**2+ &
     (vs3(1:lx1,1:lx2,1:lx3,1)-vn3)**2    !gets used several times in this subprogram



!!!!!!!!!!!!!!!!!!!!!!!!!!! O+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O+ + N2 --> NO+ + N
Teff=28/(16+28._wp)*(16*amu/3/kB*(dv2) &
     + Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn


where (Teff<=3725)
  kreac=1.71676e-12_wp &
    -7.19934e-13_wp*(Teff/300) &
    +1.33276e-13_wp*(Teff/300)**2 &
    -9.28213e-15_wp*(Teff/300)**3 &
    +6.39557e-16_wp*(Teff/300)**4
end where
where (Teff>3725 .and. Teff<=30000)
  kreac=-1.52489e-11_wp &
    +7.67112e-13_wp*(Teff/300) &
    +1.19064e-13_wp*(Teff/300)**2 &
    -1.30858e-15_wp*(Teff/300)**3 &
    +4.67756e-18_wp*(Teff/300)**4
end where
where (Teff>30000)
  kreac=-1.52489e-11_wp &
    +7.67112e-13_wp*(100) &
    +1.19064e-13_wp*(100)**2 &
    -1.30858e-15_wp*(100)**3 &
    +4.67756e-18_wp*(100)**4
end where

betanow=kreac*nn(:,:,:,2)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + O2 --> O2+ + O
Teff=32/(16+32._wp)*(16*amu/3/kB*(dv2) &
     +Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn

where (Teff<=4800)
  kreac=2.78932e-11_wp &
    -6.92612e-12_wp*(Teff/300) &
    +8.67684e-13_wp*(Teff/300)**2 &
    -3.47251e-14_wp*(Teff/300)**3 &
    +5.07097e-16_wp*(Teff/300)**4
end where
where (Teff>4800 .and. Teff<=30000)
  kreac=-1.74046e-11_wp &
    +3.02328e-12_wp*Teff/300 &
    -2.39214e-15_wp*(Teff/300)**2 &
    -4.02394e-17_wp*(Teff/300)**3
end where
where(Teff>30000)
  kreac=-1.74046e-11_wp &
    +3.02328e-12_wp*100 &
    -2.39214e-15_wp*100**2 &
    -4.02394e-17_wp*100**3
end where

betanow=kreac*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + NO --> NO+ + O
Teff=30/(16+30._wp)*(16*amu/3/kB*(dv2) &
       + Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn

where (Teff<=3800)
  kreac=6.40408e-13_wp &
    -1.33888e-13_wp*(Teff/300) &
    +7.65103e-14_wp*(Teff/300)**2 &
    -3.11509e-15_wp*(Teff/300)**3 &
    +6.62374e-17_wp*(Teff/300)**4
end where
where (Teff>3800 .and. Teff<=30000)
  kreac=-7.48312e-13_wp &
    +2.31502e-13_wp*(Teff/300) &
    +3.07160e-14_wp*(Teff/300)**2 &
    -2.65436e-16_wp*(Teff/300)**3 &
    +7.76665e-19_wp*(Teff/300)**4
end where
where (Teff>30000)
  kreac=-7.48312e-13_wp &
    +2.31502e-13_wp*(100) &
    +3.07160e-14_wp*(100)**2 &
    -2.65436e-16_wp*(100)**3 &
    +7.76665e-19_wp*(100)**4
end where

betanow=kreac*nn(:,:,:,6)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + e --> O + hv
betanow=3.7e-12_wp*(250/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.7*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!N2+ + O --> O+ + N2
Teff=16/(28+16._wp)*(28*amu/3/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) -Tn) + Tn

where (Teff <= 1500)
  kreac=1e-11_wp*(300/Teff)**0.23
elsewhere
  kreac=3.6e-12_wp*(300/Teff)**(-0.41)
end where

betanow=kreac*nn(:,:,:,1)*1e-6_wp
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N+ + O --> O+ + N
betanow=5e-13_wp*nn(:,:,:,1)*1e-6_wp
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!H+ + O --> O+ + H
Teff=Ts(1:lx1,1:lx2,1:lx3,6)
betanow = (6.e-10_wp)*(8/9._wp)*(((Teff+Tn/4)/(Tn+Teff/16))**0.5)*nn(:,:,:,1)*1e-6_wp
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,6)
Lo(:,:,:,6)=Lo(:,:,:,6)+betanow


!O+ + H --> H+ + O
betanow = 6.0e-10_wp*nn(:,:,:,4)*1e-6_wp
Pr(:,:,:,6)=Pr(:,:,:,6)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! NO+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O+ + NO --> NO+ + O Above


!O2+ + N2 --> NO+ + NO
betanow=5e-16_wp*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!O2+ + N --> NO+ + O
betanow=1.2e-10_wp*nn(:,:,:,5)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!O2+ + NO --> NO+ + O2
betanow=4.6e-10_wp*nn(:,:,:,6)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!N2+ + O --> NO+ + N
Teff=16/(28+16._wp)*(28*amu/3/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) - Tn) + Tn

where (Teff <= 1500)
  kreac=1.4e-10_wp*(300/Teff)**0.44
elsewhere
  kreac=5.2e-11_wp*(300/Teff)**(-0.2)
end where

betanow=kreac*nn(:,:,:,1)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N2+ + NO --> NO+ + N2
betanow=4.1e-10_wp*nn(:,:,:,6)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N+ + O2 --> NO+ + O
betanow=2.6e-10_wp*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!NO+ + e --> N + O
betanow=4.2e-7_wp*(300/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.85*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,2)=Lo(:,:,:,2)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! N2+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!N2+ + O2 --> O2+ + N2
Teff=32/(28+32._wp)*(28*amu/3/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) - Tn) + Tn

betanow=5e-11_wp*(300/Teff)*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N2+ + O --> NO+ + N   Above


!N2+ + O --> O+ + N2    Above


!N2+ + O --> NO+ + N    Above


!N2+ + e --> N + N
betanow=1.8e-7_wp*(300/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.39*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! O2+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O2+ + NO --> NO+ + O2  Above


!O+ + O2 --> O2+ + O    Above


!N2+ + O2 --> O2+ + N2  Above


!N+ + O2 --> O2+ + N
betanow=3.1e-10_wp*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!O2+ + e- --> O + O
where (Te <= 1200)
  kreac=1.95e-7_wp* (300/Te)**0.70! See idl code. this may need another te term
elsewhere
  kreac=7.38e-8_wp*(1200/Te)**0.56! See idl code. this may need another te term
end where

betanow=kreac*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! N+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!N+ + O --> O+ + N  Above


!N+ + O2 --> NO+ + O Above


!N+ + O2 --> O2+ + N    Above


!N+ + H --> H+ + N
betanow = 3.6e-12_wp*nn(:,:,:,4)*1e-6_wp
Pr(:,:,:,6)=Pr(:,:,:,6)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! H+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!H+ + O --> O+ + H above


!O+ + H --> H+ + O above


!N+ + H --> H+ + N above


!H+ + e --> H + hv
betanow=3.7e-12_wp*(250/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.7*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,6)=Lo(:,:,:,6)+betanow

end subroutine srcsContinuity


subroutine srcsMomentum_curv(nn,vn1,Tn,ns,vs1,vs2,vs3,Ts,E1,Q,x,Pr,Lo)

!------------------------------------------------------------
!-------POPULATE SOURCE/LOSS ARRAYS FOR MOMENTUM EQUATION.  ION
!-------PARAMETER ARGUMENTS (AND GRID STUFF) SHOULD INCLUDE GHOST CELLS
!-------NOTE THAT THIS IS THE ONLY SOURCE SUBPROGRAM WHOSE CODE
!-------DIFFERS FROM CARTESIAN TO CURVILINEAR DUE TO PRESSURE
!-------GRADIENT.
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: vn1,Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts
real(wp), dimension(-1:,-1:,-1:), intent(in) :: E1
real(wp), dimension(:,:,:,:), intent(in) :: Q
class(curvmesh), intent(in) :: x

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(inout) :: Pr,Lo
!! intent(out)

integer :: lx1,lx2,lx3,isp,isp2
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: nu,Phisj,Psisj
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: pressure,gradlp1,Epol1,gradQ
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: h1h2h3
real(wp), dimension(0:size(Ts,1)-3,size(Ts,2)-4,size(Ts,3)-4) :: tmpderiv
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: dh2dx1,dh3dx1,geom
real(wp), dimension(size(E1,1)-4,size(E1,2)-4,size(E1,3)-4) :: E1filt

integer :: ix1,ix2,ix3


lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4

Pr=0
Lo=0


!CALCULATE COMMON GEOMETRIC FACTORS USED IN EACH OF THE SPECIES CALCULATIONS
h1h2h3=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
tmpderiv=grad3D1(x%h2(0:lx1+1,1:lx2,1:lx3),x,0,lx1+1,1,lx2,1,lx3)
dh3dx1=tmpderiv(1:lx1,1:lx2,1:lx3)
tmpderiv=grad3D1(x%h3(0:lx1+1,1:lx2,1:lx3),x,0,lx1+1,1,lx2,1,lx3)
dh2dx1=tmpderiv(1:lx1,1:lx2,1:lx3)


!AMBIPOLAR ELECTRIC FIELD
pressure=ns(1:lx1,1:lx2,1:lx3,lsp)*kB*Ts(1:lx1,1:lx2,1:lx3,lsp)
gradlp1=grad3D1(log(pressure),x,1,lx1,1,lx2,1,lx3)
Epol1=kB*Ts(1:lx1,1:lx2,1:lx3,lsp)/qs(lsp)*gradlp1


!THE FIELD INTEGRATED SOLVE ELECTRIC FIELDS ARE NOT RELIABLE BELOW 100KM - AT LEAST NOT ENOUGH TO USE IN THIS CALCULATION
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if (x%alt(ix1,ix2,ix3)<100e3_wp) then
        E1filt(ix1,ix2,ix3)=0
      else
        E1filt(ix1,ix2,ix3)=E1(ix1,ix2,ix3)
      end if
    end do
  end do
end do


do isp=1,lsp
  !ION-NEUTRAL COLLISIONS
  do isp2=1,ln
    call maxwell_colln(isp,isp2,nn,Tn,Ts,nu)

    Lo(:,:,:,isp)=Lo(:,:,:,isp)+nu
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*nu*vn1
  end do


  !ION-ION
  do isp2=1,lsp
    call coulomb_colln(isp,isp2,ns,Ts,vs1,nu,Phisj,Psisj)

    Lo(:,:,:,isp)=Lo(:,:,:,isp)+nu*Phisj
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp) &
                  *nu*Phisj*vs1(1:lx1,1:lx2,1:lx3,isp2)
  end do


  !ION PRESSURE
  pressure=ns(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)
  gradlp1=grad3D1(log(pressure),x,1,lx1,1,lx2,1,lx3)                         !derivative should be from 1:lx1
  !might need to limit the gradient to non-null points like 2D MATLAB code


  !ARTIFICIAL VISCOSITY
  gradQ=grad3D1(Q(:,:,:,isp),x,1,lx1,1,lx2,1,lx3)                         !derivative should be from 1:lx1


  !GEOMETRIC FACTORS ARISING FROM ADVECTINO OF 1-COMPONENT OF MOMENTUM DENSITY
  geom=(vs2(1:lx1,1:lx2,1:lx3,isp)**2*x%h3(1:lx1,1:lx2,1:lx3)*dh2dx1+ &
        vs3(1:lx1,1:lx2,1:lx3,isp)**2*x%h2(1:lx1,1:lx2,1:lx3)*dh3dx1)*ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)/h1h2h3


  !ACCUMULATED ALL FORCES
!      Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)*(E1+Epol1) &
  Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)*(E1filt+Epol1) &
                -pressure*gradlp1 &
                -gradQ &
                +geom &
                +ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*x%g1
end do

end subroutine srcsMomentum_curv


subroutine srcsEnergy(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts,Pr,Lo,E1,E2,E3,x)

!------------------------------------------------------------
!-------POPULATE SOURCE/LOSS ARRAYS FOR ENERGY EQUATION.  ION
!-------PARAMETER ARGUMENTS SHOULD INCLUDE GHOST CELLS
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts
real(wp), dimension(-1:,-1:,-1:), intent(in) :: E1,E2,E3
class(curvmesh), intent(in) :: x !Added for FBI, need BMAG


real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(inout) :: Pr,Lo
!! intent(out)

integer :: ix1,ix2,ix3,lx1,lx2,lx3,isp,isp2
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: nu,Phisj,Psisj
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: fact,iePT,ieLT,f,g    !work arrays
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: iePTFBI,ieLTFBI    !FBI array
real(wp) :: sfact

lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4

Pr=0
Lo=0
iePT=0
ieLT=0


!ELASTIC COLLISIONS
do isp=1,lsp
  !ION-NEUTRAL
  do isp2=1,ln
    call maxwell_colln(isp,isp2,nn,Tn,Ts,nu)    ! isp=7 is electrons

    !HEAT TRANSFER
    fact=2*nu/(ms(isp)+mn(isp2))
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*kB/(gammas(isp)-1)*fact*Tn
    Lo(:,:,:,isp)=Lo(:,:,:,isp)+ms(isp)*fact


    !FRICTION
    fact=fact*mn(isp2)/3
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)/(gammas(isp)-1) &
                  *((vs1(1:lx1,1:lx2,1:lx3,isp)-vn1)**2+(vs2(1:lx1,1:lx2,1:lx3,isp)-vn2)**2 &
                  +(vs3(1:lx1,1:lx2,1:lx3,isp)-vn3)**2)*fact     !vn's should be correct shape for this...
  end do

  !ION-ION
  do isp2=1,lsp
    call coulomb_colln(isp,isp2,ns,Ts,vs1,nu,Phisj,Psisj)

    !HEAT TRANSFER
    fact=2*nu*Psisj/(ms(isp)+ms(isp2))
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*kB/(gammas(isp)-1) &
                  *fact*Ts(1:lx1,1:lx2,1:lx3,isp2)
    Lo(:,:,:,isp)=Lo(:,:,:,isp)+ms(isp)*fact

    !FRICTION
!        fact=2*nu*Phisj/(ms(isp)+ms(isp2))*mn(isp2)/3     !this is the error that was causing the runtime problem with -O3 on phys_consts.f90.  Much thanks to Guy Grubbs for finding this longstanding error.
    fact=2*nu*Phisj/(ms(isp)+ms(isp2))*ms(isp2)/3
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)/(gammas(isp)-1) &
                  *((vs1(1:lx1,1:lx2,1:lx3,isp)-vs1(1:lx1,1:lx2,1:lx3,isp2))**2 &
                   +(vs2(1:lx1,1:lx2,1:lx3,isp)-vs2(1:lx1,1:lx2,1:lx3,isp2))**2 &
                   +(vs3(1:lx1,1:lx2,1:lx3,isp)-vs3(1:lx1,1:lx2,1:lx3,isp2))**2)*fact
  end do
end do


!INELASTIC COLLISIONS FOR ELECTRONS, ROTATIONAL
sfact=elchrg/kB*(gammas(lsp)-1);   !cf. S&N 2010, electron energy equatoin section
nu=sfact*6.9e-14_wp*nn(:,:,:,3)*1e-6_wp/sqrt(Ts(1:lx1,1:lx2,1:lx3,lsp))    !O2 rotational excitation
iePT=nu*Tn
ieLT=nu
nu=sfact*2.9e-14_wp*nn(:,:,:,2)*1e-6_wp/sqrt(Ts(1:lx1,1:lx2,1:lx3,lsp))    !N2 rot. exc.
iePT=iePT+nu*Tn;
ieLT=ieLT+nu;

f=1.06e4_wp+7.51e3_wp*tanh(1.10e-3_wp*(Ts(1:lx1,1:lx2,1:lx3,lsp)-1800))
g=3300+1.233_wp*(Ts(1:lx1,1:lx2,1:lx3,lsp)-1000)-2.056e-4_wp &
  *(Ts(1:lx1,1:lx2,1:lx3,lsp)-1000)*(Ts(1:lx1,1:lx2,1:lx3,lsp)-4000)
fact=sfact*2.99e-12_wp*nn(:,:,:,2)*1e-6_wp*exp(f*(Ts(1:lx1,1:lx2,1:lx3,lsp)-2000) &
        /Ts(1:lx1,1:lx2,1:lx3,lsp)/2000)*(exp(-g*(Ts(1:lx1,1:lx2,1:lx3,lsp)-Tn) &
        /Ts(1:lx1,1:lx2,1:lx3,lsp)/Tn)-1)    !N2 vibrational excitation
iePT=iePT-max(fact,0._wp);
f=3300-839*sin(1.91e-4_wp*(Ts(1:lx1,1:lx2,1:lx3,lsp)-2700))
fact=sfact*5.196e-13_wp*nn(:,:,:,3)*1e-6_wp*exp(f*(Ts(1:lx1,1:lx2,1:lx3,lsp)-700) &
     /Ts(1:lx1,1:lx2,1:lx3,lsp)/700)*(exp(-2770*(Ts(1:lx1,1:lx2,1:lx3,lsp)-Tn) &
     /Ts(1:lx1,1:lx2,1:lx3,lsp)/Tn)-1)    !O2 vibrational excitation
iePT=iePT-max(fact,0._wp);


!! This would be the place to include FBI heating probably just add to iePT
! iePT=iePT+FBIheating()
call FBIheating(nn,Tn,ns,Ts,E1,E2,E3,x,iePTFBI,IeLTFBI)
iePT=iePT+iePTFBI
ieLT=ieLT*ieLTFBI 
!now, here is have questions, since I do not how really what this terms does (loss factor), where to apply it, or if this is the right place


!CORRECT TEMP EXPRESSIONS TO CORRESPOND TO INTERNAL ENERGY SOURCES
Pr(:,:,:,lsp)=Pr(:,:,:,lsp)+iePT*ns(1:lx1,1:lx2,1:lx3,lsp)*kB/(gammas(lsp)-1)   !Arg, forgot about the damn ghost cells in original code...
Lo(:,:,:,lsp)=Lo(:,:,:,lsp)+ieLT

end subroutine srcsEnergy




subroutine FBIheating(nn,Tn,ns,Ts,E1,E2,E3,x,iePTFBI,ieLTFBI)

!! Inputs Needed
real(wp), dimension(:,:,:,:), intent(in) :: nn !Neutral density
real(wp), dimension(:,:,:), intent(in) :: Tn !neutral temperature
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts !Plasma density and temperature
real(wp), dimension(-1:,-1:,-1:), intent(in) :: E1,E2,E3 !Electric Field
class(curvmesh), intent(in) :: x !Grid, doing this because BMAG is stored here

!! intent(out)
real(wp), dimension(:,:,:), intent(inout) :: iePTFBI,ieLTFBI ! Two terms, one is heating and the other one is a factor for cooling. 

!!Internal Arrays
integer :: isp,isp2,lx1,lx2,lx3
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp) :: nsuAvg 
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp-1) :: niW 
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,ln) :: nuW 
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,2) :: nuAvg, msAvg, TsAvg  
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: Bmagnitude, nu, nsAvg, omegae, omegai, ki, ke, phi
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: Ethresholdnum, Ethresholdden, Ethreshold, Emagnitude
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: heatingfirst, heatingsecond, heatingtotal, lossfactor

!I have a lot of internal variables, some might be redundant
!nsuAvg= Collision frequencies for each species averaged over all neutrals
!nuW+ MassDensity Weight for averaging
!nuAvg= Collision frequencies for ions (:,:,:,1) averaged over al species and electrons (:,:,:,2)
!msAvg= Average ion mass (:,:,:,1), wieghted by the density at each pixel. Electron mass is also included (:,:,:,2). No Ghost cells
!TsAvg= Average ion temperature (:,:,:,1), wieghted by the density at each pixel. Electron temperature is also included (:,:,:,2) No Ghost cells
!Bmagnitude= Magnitude of the magnetic field WITHOUT ghost cells
!nu= Output of maxwell_coll, used inside a loop
!nsAvg= Average density of each pixel, currently it assumes qneutrality, could be wrong. 
!omegai= Plasma parameter omega for ions
!omegai= Plasma parameter omega for electrons
!ki and ke=Omega diveded by collision frequency, used in FBI heating equation
!phi= Inverse of the product between ki and ke, used in FBI heating equation
!Ethresholdnum= Numeretor of the Ethreshold equation, done this way to avoid mistakes
!Ethresholden= Denominator of the Ethreshold equation, done this way to avoid mistakes
!Ethreshold= Threshold at which FBI should appear (this might be a redundant thing, could eliminate the two from above if works correctly)
!Emagnitude= Magnitude of the electric field WITHOUT ghost cells
!heatingfirst= First term of the FBI heating equation
!heatingsecond= Second term of the FBI heating equation
!heatingtotal= Multiplication of both terms from above, again can be made into just one variuable if this works.
!Lossfactor= Factor used to reduce supertherman electron velocities, ergo reduced electron cooling rate.  
lx1=x%lx1
lx2=x%lx2
lx3=x%lx3

!Bmagnitude=x%Bmag(1:lx1,1:lx2,1:lx3)
Bmagnitude=x%Bmag(1:lx1,1:lx2,1:lx3)
print*,shape(E1),shape(E2),shape(E3)
Emagnitude=sqrt(E1(1:lx1,1:lx2,1:lx3)**2+E2(1:lx1,1:lx2,1:lx3)**2+E3(1:lx1,1:lx2,1:lx3)**2) !!Already evaluated with no ghost cells
!!Indices for stuff
lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4

!!Initialize arrays a 0s and loss as 1s
nuAvg=0.0
nsuAvg=0.0
msAvg=0.0
nsAvg=0.0
TsAvg=0.0
iePTFBI=0.0
ieLTFBI=0.0
lossfactor=1.0

!! Find average collision frequencies. Loop over all neutrals, store ion and electric collision frequencies. Use average?
!! Just store all of them, decide later. 

!MassDensity Weight of Neutrals
do isp2=1,ln
    nuW(:,:,:,isp2)=nn(:,:,:,isp2)*mn(isp2) !Wieght of the neutrals
end do

!!MassDensity Weight of ions
do isp=1,lsp-1
  niW(:,:,:,isp)=ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)
end do

!! Average Collusuons frequencies: loop for doing average of collisions frequencies, first averaging over neutrals
do isp=1,lsp
  do isp2=1,ln
    call maxwell_colln(isp,isp2,nn,Tn,Ts,nu)
    nsuAvg(:,:,:,isp)=nsuAvg(:,:,:,isp)+nu*nuW(:,:,:,isp2) !Store the collision frequencies weighted by massdensity
  end do
  nsuAvg(:,:,:,isp)=nsuAvg(:,:,:,isp)/sum(nuW, dim=4) !! Average over all neutrals weighted by MassDensity
end do

!!Loop for average over ions
do isp=1,lsp-1
  nuAvg(:,:,:,1)=nuAvg(:,:,:,1)+nsuAvg(:,:,:,isp)*niW(:,:,:,isp) !! Store summation of collision frequencies weighted by MassDensity
end do
nuAvg(:,:,:,1)=nuAvg(:,:,:,1)/sum(niW,dim=4) !! Average over all ions weighted by MassDensity

!!Electrons do not need averaging
nuAvg(:,:,:,2)=nsuAvg(:,:,:,lsp)

!! Average mass of ions, also weighted by MassDensity
do isp=1,lsp-1
  msAvg(:,:,:,1)=msAvg(:,:,:,1)+ms(isp)*niW(:,:,:,isp)
end do

msAvg(:,:,:,1)=msAvg(:,:,:,1)/(lsp-1)/sum(niW,dim=4)

msAvg(:,:,:,2)=ms(lsp) !! Electron mass

!! Average density
nsAvg=ns(1:lx1,1:lx2,1:lx3,lsp) !! Assume qneutrality, could be wrong

!! ki value
omegai=elchrg*Bmagnitude/msAvg(:,:,:,1) !! Would this work?, it will, I defined Bmagnitude above
ki=abs(omegai/nuAvg(:,:,:,1)) !! Could do ABS to be sure of the sign
!! ke value
omegae=elchrg*Bmagnitude/msAvg(:,:,:,2) 
ke=abs(omegae/nuAvg(:,:,:,2))

!!Phi value 1/(ki*ke)
phi=1/(ke*ki)

!!Average ion temperature
do isp=1,lsp-1
  TsAvg(:,:,:,1)=TsAvg(:,:,:,1)+Ts(1:lx1,1:lx2,1:lx3,isp)*niW(:,:,:,isp)
end do
TsAvg(:,:,:,1)=TsAvg(:,:,:,1)/sum(niW,dim=4) 
TsAvg(:,:,:,2)=Ts(1:lx1,1:lx2,1:lx3,lsp)

!!Ethreshold
Ethresholdnum=(1+phi)*SQRT(kB*(1+ki**2)*(TsAvg(:,:,:,1)+TsAvg(:,:,:,2))*Bmagnitude)
Ethresholdden=sqrt((1-ki**2)*msAvg(:,:,:,1)) 

Ethreshold=Ethresholdnum(:,:,:)/Ethresholdden(:,:,:)

!!First term of heating equation, have to check the Emagnitude part
heatingfirst=msAvg(:,:,:,1)*nuAvg(:,:,:,1)*nsAvg*(ki**2)*(Emagnitude-Ethreshold)**2/((1+ki**2)*Bmagnitude**2)
heatingsecond=(Emagnitude*(1+phi)/Ethreshold-1)

!heating total as if FBI was everywhere
heatingtotal=heatingfirst*heatingsecond

!Loss factor a every pixel, as if FBI was everywhere
lossfactor=exp(-7.54e-4_wp*(TsAvg(:,:,:,2)-500))

where (Emagnitude<Ethreshold) !Anything without a sufficiente E field gets back to normal.
  heatingtotal=0.0
  lossfactor=1.0
end where

where (ki>=1) !Anything where ions are magnetized also goes back to normal
  heatingtotal=0.0
  lossfactor=1.0
end where

iePTFBI=heatingtotal
ieLTFBI=lossfactor

end subroutine FBIheating


end module sources

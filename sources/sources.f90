module sources

use phys_consts, only: wp, lsp, amu, kb, qs, ln, ms, gammas, elchrg, mn
use grid, only : curvmesh,g1
use calculus, only : grad3d1
use mpimod, only: myid, tagvs1bc, tagvs2bc, tagvs3bc, lid, halo23
use collisions, only:  maxwell_colln, coulomb_colln

implicit none

interface srcsMomentum
  module procedure srcsMomentum_curv
end interface srcsMomentum

private
public ::  rk2_prep_mpi, srcsenergy, srcsmomentum, srcscontinuity

contains


subroutine srcsContinuity(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts,Pr,Lo)

!------------------------------------------------------------
!-------POPULATE SOURCE/LOSS ARRAYS FOR CONTINUITY EQUATION.  ION
!-------PARAMETER ARGUMENTS (AND GRID STUFF) SHOULD INCLUDE GHOST CELLS
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(out) :: Pr,Lo

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: betanow,kreac,Teff,Te,dv2
integer :: lx1,lx2,lx3

lx1=size(ns,1)-4
lx2=size(ns,2)-4
lx3=size(ns,3)-4

Pr=0._wp
Lo=0._wp
Te=Ts(1:lx1,1:lx2,1:lx3,lsp)    !used at all???  needs to be checked and removed, if not
dv2=(vs1(1:lx1,1:lx2,1:lx3,1)-vn1)**2+(vs2(1:lx1,1:lx2,1:lx3,1)-vn2)**2+ &
     (vs3(1:lx1,1:lx2,1:lx3,1)-vn3)**2    !gets used several times in this subprogram



!!!!!!!!!!!!!!!!!!!!!!!!!!! O+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O+ + N2 --> NO+ + N
Teff=28d0/(16d0+28d0)*(16d0*amu/3d0/kB*(dv2) &
     + Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn


where (Teff<=3725d0)
  kreac=1.71676d-12 &
    -7.19934d-13*(Teff/300._wp) &
    +1.33276d-13*(Teff/300._wp)**2 &
    -9.28213d-15*(Teff/300._wp)**3 &
    +6.39557d-16*(Teff/300._wp)**4
end where
where (Teff>3725d0 .and. Teff<=30000._wp)
  kreac=-1.52489d-11 &
    +7.67112d-13*(Teff/300._wp) &
    +1.19064d-13*(Teff/300._wp)**2 &
    -1.30858d-15*(Teff/300._wp)**3 &
    +4.67756d-18*(Teff/300._wp)**4
end where
where (Teff>30000._wp)
  kreac=-1.52489d-11 &
    +7.67112d-13*(100._wp) &
    +1.19064d-13*(100._wp)**2 &
    -1.30858d-15*(100._wp)**3 &
    +4.67756d-18*(100._wp)**4
end where

betanow=kreac*nn(:,:,:,2)*1d-6
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + O2 --> O2+ + O
Teff=32._wp/(16d0+32._wp)*(16d0*amu/3d0/kB*(dv2) &
     +Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn

where (Teff<=4800._wp)
  kreac=2.78932d-11 &
    -6.92612d-12*(Teff/300._wp) &
    +8.67684d-13*(Teff/300._wp)**2 &
    -3.47251d-14*(Teff/300._wp)**3 &
    +5.07097d-16*(Teff/300._wp)**4
end where
where (Teff>4800._wp .and. Teff<=30000._wp)
  kreac=-1.74046d-11 &
    +3.02328d-12*(Teff)/300._wp &
    -2.39214d-15*(Teff/300._wp)**2 &
    -4.02394d-17*(Teff/300._wp)**3
end where
where(Teff>30000._wp)
  kreac=-1.74046d-11 &
    +3.02328d-12*(100._wp) &
    -2.39214d-15*(100._wp)**2 &
    -4.02394d-17*(100._wp)**3
end where

betanow=kreac*nn(:,:,:,3)*1d-6
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + NO --> NO+ + O
Teff=30._wp/(16d0+30._wp)*(16d0*amu/3d0/kB*(dv2) &
       + Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn

where (Teff<=3800._wp)
  kreac=6.40408d-13 &
    -1.33888d-13*(Teff/300._wp) &
    +7.65103d-14*(Teff/300._wp)**2 &
    -3.11509d-15*(Teff/300._wp)**3 &
    +6.62374d-17*(Teff/300._wp)**4
end where
where (Teff>3800._wp .and. Teff<=30000._wp)
  kreac=-7.48312d-13 &
    +2.31502d-13*(Teff/300._wp) &
    +3.07160d-14*(Teff/300._wp)**2 &
    -2.65436d-16*(Teff/300._wp)**3 &
    +7.76665d-19*(Teff/300._wp)**4
end where
where (Teff>30000._wp)
  kreac=-7.48312d-13 &
    +2.31502d-13*(100._wp) &
    +3.07160d-14*(100._wp)**2 &
    -2.65436d-16*(100._wp)**3 &
    +7.76665d-19*(100._wp)**4
end where

betanow=kreac*nn(:,:,:,6)*1d-6
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + e --> O + hv
betanow=3.7d-12*(250._wp/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.7*ns(1:lx1,1:lx2,1:lx3,lsp)*1d-6
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!N2+ + O --> O+ + N2
Teff=16d0/(28d0+16d0)*(28d0*amu/3d0/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) -Tn) + Tn

where (Teff <= 1500._wp)
  kreac=1d-11*(300._wp/Teff)**0.23
elsewhere
  kreac=3.6d-12*(300._wp/Teff)**(-0.41)
end where

betanow=kreac*nn(:,:,:,1)*1d-6
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N+ + O --> O+ + N
betanow=5d-13*nn(:,:,:,1)*1d-6
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!H+ + O --> O+ + H
Teff=Ts(1:lx1,1:lx2,1:lx3,6)
betanow = (6.0d-10)*(8d0/9d0)*(((Teff+Tn/4d0)/(Tn+Teff/16d0))**0.5)*nn(:,:,:,1)*1d-6
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,6)
Lo(:,:,:,6)=Lo(:,:,:,6)+betanow


!O+ + H --> H+ + O
betanow = 6.0d-10*nn(:,:,:,4)*1d-6
Pr(:,:,:,6)=Pr(:,:,:,6)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! NO+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O+ + NO --> NO+ + O Above


!O2+ + N2 --> NO+ + NO
betanow=5d-16*nn(:,:,:,3)*1d-6
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!O2+ + N --> NO+ + O
betanow=1.2d-10*nn(:,:,:,5)*1d-6
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!O2+ + NO --> NO+ + O2
betanow=4.6d-10*nn(:,:,:,6)*1d-6
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!N2+ + O --> NO+ + N
Teff=16d0/(28d0+16d0)*(28d0*amu/3d0/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) - Tn) + Tn

where (Teff <= 1500._wp)
  kreac=1.4d-10*(300._wp/Teff)**0.44
elsewhere
  kreac=5.2d-11*(300._wp/Teff)**(-0.2)
end where

betanow=kreac*nn(:,:,:,1)*1d-6
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N2+ + NO --> NO+ + N2
betanow=4.1d-10*nn(:,:,:,6)*1d-6
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N+ + O2 --> NO+ + O 
betanow=2.6d-10*nn(:,:,:,3)*1d-6
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!NO+ + e --> N + O 
betanow=4.2d-7*(300._wp/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.85*ns(1:lx1,1:lx2,1:lx3,lsp)*1d-6
Lo(:,:,:,2)=Lo(:,:,:,2)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! N2+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!N2+ + O2 --> O2+ + N2
Teff=32._wp/(28d0+32._wp)*(28d0*amu/3d0/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) - Tn) + Tn

betanow=5d-11*(300._wp/Teff)*nn(:,:,:,3)*1d-6
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N2+ + O --> NO+ + N   Above


!N2+ + O --> O+ + N2    Above


!N2+ + O --> NO+ + N    Above


!N2+ + e --> N + N
betanow=1.8d-7*(300._wp/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.39*ns(1:lx1,1:lx2,1:lx3,lsp)*1d-6
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! O2+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O2+ + NO --> NO+ + O2  Above


!O+ + O2 --> O2+ + O    Above


!N2+ + O2 --> O2+ + N2  Above


!N+ + O2 --> O2+ + N
betanow=3.1d-10*nn(:,:,:,3)*1d-6
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!O2+ + e- --> O + O
where (Te <= 1200._wp)
  kreac=1.95d-7* (300._wp/Te)**0.70! See idl code. this may need annother te term
elsewhere
  kreac=7.38d-8*(1200/Te)**0.56! See idl code. this may need annother te term
end where

betanow=kreac*ns(1:lx1,1:lx2,1:lx3,lsp)*1d-6
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! N+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!N+ + O --> O+ + N  Above


!N+ + O2 --> NO+ + O Above


!N+ + O2 --> O2+ + N    Above


!N+ + H --> H+ + N
betanow = 3.6d-12*nn(:,:,:,4)*1d-6
Pr(:,:,:,6)=Pr(:,:,:,6)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! H+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!H+ + O --> O+ + H above


!O+ + H --> H+ + O above


!N+ + H --> H+ + N above


!H+ + e --> H + hv
betanow=3.7d-12*(250/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.7*ns(1:lx1,1:lx2,1:lx3,lsp)*1d-6
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
real(wp), dimension(:,:,:), intent(in) :: E1
real(wp), dimension(:,:,:,:), intent(in) :: Q
type(curvmesh), intent(in) :: x

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(out) :: Pr,Lo

integer :: lx1,lx2,lx3,isp,isp2
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: nu,Phisj,Psisj
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: pressure,gradlp1,Epol1,gradQ
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: h1h2h3
real(wp), dimension(0:size(Ts,1)-3,size(Ts,2)-4,size(Ts,3)-4) :: tmpderiv
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: dh2dx1,dh3dx1,geom
real(wp), dimension(size(E1,1),size(E1,2),size(E1,3)) :: E1filt

integer :: ix1,ix2,ix3


lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4

Pr=0._wp
Lo=0._wp


!CALCULATE COMMON GEOMETRIC FACTORS USED IN EACH FO THE SPECIES CALCULATIONS
h1h2h3=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
tmpderiv=grad3D1(x%h2(0:lx1+1,1:lx2,1:lx3),x,0,lx1+1,1,lx2,1,lx3)
dh3dx1=tmpderiv(1:lx1,1:lx2,1:lx3)
tmpderiv=grad3D1(x%h3(0:lx1+1,1:lx2,1:lx3),x,0,lx1+1,1,lx2,1,lx3)
dh2dx1=tmpderiv(1:lx1,1:lx2,1:lx3)


!AMBIPOLAR ELECTRIC FIELD
pressure=ns(1:lx1,1:lx2,1:lx3,lsp)*kB*Ts(1:lx1,1:lx2,1:lx3,lsp)
gradlp1=grad3D1(log(pressure),x,1,lx1,1,lx2,1,lx3)
Epol1=kB*Ts(1:lx1,1:lx2,1:lx3,lsp)/qs(lsp)*gradlp1


!THE FIELD INTEGRATED SOLVE ELECTRIC FIELDS ARE NOT RELIABLE BELOW 100KM - AT LEAST NOT ENOGUH TO USE IN THIS CALCULATION
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if (x%alt(ix1,ix2,ix3)<100d3) then
        E1filt(ix1,ix2,ix3)=0._wp
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
                +ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*g1
end do

end subroutine srcsMomentum_curv


subroutine srcsEnergy(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts,Pr,Lo)

!------------------------------------------------------------
!-------POPULATE SOURCE/LOSS ARRAYS FOR ENERGY EQUATION.  ION
!-------PARAMETER ARGUMENTS SHOULD INCLUDE GHOST CELLS
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(out) :: Pr,Lo

integer :: ix1,ix2,ix3,lx1,lx2,lx3,isp,isp2
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: nu,Phisj,Psisj
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: fact,iePT,ieLT,f,g    !work arrays
real(wp) :: sfact

lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4

Pr=0._wp
Lo=0._wp
iePT=0._wp
ieLT=0._wp


!ELASTIC COLLISIONS
do isp=1,lsp
  !ION-NEUTRAL
  do isp2=1,ln
    call maxwell_colln(isp,isp2,nn,Tn,Ts,nu)

    !HEAT TRANSFER
    fact=2._wp*nu/(ms(isp)+mn(isp2))
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*kB/(gammas(isp)-1._wp)*fact*Tn
    Lo(:,:,:,isp)=Lo(:,:,:,isp)+ms(isp)*fact
   

    !FRICTION (neglect vn for now)
    fact=fact*mn(isp2)/3d0
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)/(gammas(isp)-1._wp) &
                  *((vs1(1:lx1,1:lx2,1:lx3,isp)-vn1)**2+(vs2(1:lx1,1:lx2,1:lx3,isp)-vn2)**2 &
                  +(vs3(1:lx1,1:lx2,1:lx3,isp)-vn3)**2)*fact     !vn's should be correct shape for this...
  end do

  !ION-ION
  do isp2=1,lsp
    call coulomb_colln(isp,isp2,ns,Ts,vs1,nu,Phisj,Psisj)

    !HEAT TRANSFER
    fact=2._wp*nu*Psisj/(ms(isp)+ms(isp2))
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*kB/(gammas(isp)-1._wp) &
                  *fact*Ts(1:lx1,1:lx2,1:lx3,isp2)
    Lo(:,:,:,isp)=Lo(:,:,:,isp)+ms(isp)*fact

    !FRICTION
!        fact=2._wp*nu*Phisj/(ms(isp)+ms(isp2))*mn(isp2)/3d0     !this is the error that was causing the runtime problem with -O3 on phys_consts.f90.  Much thanks to Guy Grubbs for finding this longstanding error.
    fact=2._wp*nu*Phisj/(ms(isp)+ms(isp2))*ms(isp2)/3d0
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)/(gammas(isp)-1._wp) &
                  *((vs1(1:lx1,1:lx2,1:lx3,isp)-vs1(1:lx1,1:lx2,1:lx3,isp2))**2 &
                   +(vs2(1:lx1,1:lx2,1:lx3,isp)-vs2(1:lx1,1:lx2,1:lx3,isp2))**2 &
                   +(vs3(1:lx1,1:lx2,1:lx3,isp)-vs3(1:lx1,1:lx2,1:lx3,isp2))**2)*fact
  end do
end do


!INELASTIC COLLISIONS FOR ELECTRONS, ROTATIONAL
sfact=elchrg/kB*(gammas(lsp)-1._wp);   !cf. S&N 2010, electron energy equatoin section
nu=sfact*6.9d-14*nn(:,:,:,3)*1d-6/sqrt(Ts(1:lx1,1:lx2,1:lx3,lsp))    !O2 rotational excitation
iePT=nu*Tn
ieLT=nu
nu=sfact*2.9d-14*nn(:,:,:,2)*1d-6/sqrt(Ts(1:lx1,1:lx2,1:lx3,lsp))    !N2 rot. exc.
iePT=iePT+nu*Tn;
ieLT=ieLT+nu;

f=1.06d4+7.51d3*tanh(1.10d-3*(Ts(1:lx1,1:lx2,1:lx3,lsp)-1800._wp))
g=3300._wp+1.233d0*(Ts(1:lx1,1:lx2,1:lx3,lsp)-1000._wp)-2.056d-4 &
  *(Ts(1:lx1,1:lx2,1:lx3,lsp)-1000._wp)*(Ts(1:lx1,1:lx2,1:lx3,lsp)-4000._wp)
fact=sfact*2.99d-12*nn(:,:,:,2)*1d-6*exp(f*(Ts(1:lx1,1:lx2,1:lx3,lsp)-2000._wp) &
        /Ts(1:lx1,1:lx2,1:lx3,lsp)/2000._wp)*(exp(-1._wp*g*(Ts(1:lx1,1:lx2,1:lx3,lsp)-Tn) &
        /Ts(1:lx1,1:lx2,1:lx3,lsp)/Tn)-1._wp)    !N2 vibrational excitation
iePT=iePT-max(fact,0._wp);
f=3300._wp-839d0*sin(1.91d-4*(Ts(1:lx1,1:lx2,1:lx3,lsp)-2700._wp))
fact=sfact*5.196d-13*nn(:,:,:,3)*1d-6*exp(f*(Ts(1:lx1,1:lx2,1:lx3,lsp)-700._wp) &
     /Ts(1:lx1,1:lx2,1:lx3,lsp)/700._wp)*(exp(-2770._wp*(Ts(1:lx1,1:lx2,1:lx3,lsp)-Tn) &
     /Ts(1:lx1,1:lx2,1:lx3,lsp)/Tn)-1._wp)    !O2 vibrational excitation
iePT=iePT-max(fact,0._wp);


!CORRECT TEMP EXPRESSIONS TO CORRESPOND TO INTERNAL ENERGY SOURCES
Pr(:,:,:,lsp)=Pr(:,:,:,lsp)+iePT*ns(1:lx1,1:lx2,1:lx3,lsp)*kB/(gammas(lsp)-1._wp)   !Arg, forgot about the damn ghost cells in original code...
Lo(:,:,:,lsp)=Lo(:,:,:,lsp)+ieLT

end subroutine srcsEnergy


subroutine RK2_prep_mpi(isp,flagperiodic,vs1,vs2,vs3)

!------------------------------------------------------------
!-------PASS BOUNDARY CELLS FOR COMPUTING COMPRESSION.  DONE ON
!-------A PER-SPECIES BASIS.  ION
!-------PARAMETER ARGUMENTS SHOULD INCLUDE GHOST CELLS.  DO
!-------WE NEED TO PASS V1,2 VARIABLES FOR DIV?
!------------------------------------------------------------

integer, intent(in) :: isp
integer, intent(in) :: flagperiodic
real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs1,vs2,vs3

real(wp), dimension(-1:size(vs1,1)-2,-1:size(vs1,2)-2,-1:size(vs1,3)-2) :: param
integer :: lx1,lx2,lx3
integer :: idleft,idright

lx1=size(vs1,1)-4
lx2=size(vs1,2)-4
lx3=size(vs1,3)-4


!ZOH EXTRAPOLATION OF V1,2 VARIABLES
vs1(0,:,:,isp)=vs1(1,:,:,isp)
vs1(lx1+1,:,:,isp)=vs1(lx1,:,:,isp)
vs2(:,0,:,isp)=vs2(:,1,:,isp)
vs2(:,lx2+1,:,isp)=vs2(:,lx2,:,isp)


!IDENTIFY MY NEIGHBORS
idleft=myid-1; idright=myid+1


!-- Now halo the interior parts (must happen for every worker since even a worker with a
!-- global boundary will still have one interior boundary to be haloed.
!BY DEFAULT THE GLOBAL BOUNDARIES ARE ASSUMED TO BE PERIOIDIC
param=vs1(:,:,:,isp)
call halo23(param,1,tagvs1BC)
vs1(:,:,:,isp)=param
param=vs2(:,:,:,isp)
call halo23(param,1,tagvs2BC)
vs2(:,:,:,isp)=param
param=vs3(:,:,:,isp)
call halo23(param,1,tagvs3BC)
vs3(:,:,:,isp)=param


!ZERO ORDER HOLD EXTRAPOLATION OF BOUNDARIES - OTHERWISE LEAVE AS PERIODIC IF REQUESTED
if (flagperiodic==0) then
  if (idleft==-1) then    !left x3 boundary
    vs1(:,:,0,isp)=vs1(:,:,1,isp)
    vs2(:,:,0,isp)=vs2(:,:,1,isp)
    vs3(:,:,0,isp)=vs3(:,:,1,isp)
  elseif (idright==lid) then    !right x3 boundary
    vs1(:,:,lx3+1,isp)=vs1(:,:,lx3,isp)
    vs2(:,:,lx3+1,isp)=vs2(:,:,lx3,isp)
    vs3(:,:,lx3+1,isp)=vs3(:,:,lx3,isp)
  end if
end if

end subroutine RK2_prep_mpi

end module sources

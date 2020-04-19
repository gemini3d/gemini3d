program MHD1D_SAW

!A 'CLASSICAL HYDRODYNAMIC' SOLVER (WAVES PROPAGATED THRU SOURCE TERMS)
!FOR THE IDEAL MHD EQUATIONS IN 1.66 DIMENSIONS.  M. ZETTERGREN, T. SYMONS

use advec
implicit none (type, external)

integer, parameter :: method=1
integer, parameter :: npts=2001                                         !hard-coded for SAW problem
real(wp), parameter :: tcfl=0.75
real(wp), parameter :: pi=3.14
real(wp), parameter :: amu=1.67e-27, kb=1.38e-23, gammap=5.0/3.0, mu0=4.0*pi*1e-7
real(wp), parameter :: stride=10e3

real(wp), dimension(npts) :: dx1i
real(wp), dimension(-1:npts+2) :: x1
real(wp), dimension(-1:npts+2) :: rhom,rhou1,rhou2,rhou3,B1,B2,B3,rhoeps,u1,u2,u3
real(wp), dimension(npts+1) :: x1i,v1i
real(wp), dimension(0:npts+2) :: dx1
integer :: lx1,it,ix1, u
real(wp) :: t=0,dt=1e-6,dtout=5,tdur=300,xi=2
real(wp) :: toutnext

real(wp), dimension(npts) :: Q,p,grad1B1u3,rhoepshalf,sourceterm,du1full,grad1u1,vA,vsnd   !'work' arrays
real(wp) :: sinwt,sinwt2                                                                   !'work' vars
real(wp) :: u31,B10,rhom0,vA0,B31,omega0,T0,tfwhm,tdel,Ti=8000                             !linear wave params



!GRID CONSTRUCTION
x1=[ (ix1*stride, ix1=-2,npts+1) ]
lx1=size(x1)-4                             !exclude ghost cells in count
dx1=x1(0:lx1+2)-x1(-1:lx1+1)                !backward diffs
x1i(1:lx1+1)=0.5*(x1(0:lx1)+x1(1:lx1+1))    !interface data
dx1i=x1i(2:lx1+1)-x1i(1:lx1)


!PREP OUTPUT FILE

open(newunit=u,file='MHD1D.dat',status='replace', access='stream', action='write')
write(u,*) lx1
call writearray(u,x1)


!INITIAL CONDITIONS
rhom=1e9*amu
rhou1=0
rhou2=0
rhou3=0
B1=100e-9
B2=0
B3=0
rhoeps=rhom*kb*Ti/amu/(gammap-1)

u1=rhou1/rhom
u2=rhou2/rhom
u3=rhou3/rhom


!STATIC INFO FOR CALCULATING BOUNDARY CONDITIONS
u31=100
B10=sum(B1)/size(B1)
rhom0=sum(rhom)/size(rhom)
vA0=B10/sqrt(mu0*rhom0)
B31=-u31*B10/vA0
omega0=0.15*2*pi
T0=2*pi/omega0
tfwhm=5*T0
tdel=10*T0


!MAIN LOOP
toutnext=dtout
do while (t<tdur)
  !ART. VISCOSITY DIFFS OF U1
  du1full=rhou1(2:lx1+1)/rhom(2:lx1+1)-rhou1(0:lx1-1)/rhom(0:lx1-1)


  !TIME STEP DETERMINATION
  vA=sqrt((B1(1:lx1)**2+B2(1:lx1)**2+B3(1:lx1)**2)/mu0/rhom(1:lx1))
  p=rhoeps(1:lx1)*(gammap-1)
  vsnd=sqrt(gammap*p/rhom(1:lx1))
  dt=tcfl*minval(dx1i/(vA+vsnd+abs(u1(1:lx1))))
  t=t+dt                                        !time after this step
  write(*,*) 'Time:  ',t,'  Time step:  ',dt, &
               '  Alfven speed:  ',maxval(vA/1e3), &
               '  sound speed:  ',maxval(vsnd/1e3), &
               '  matter speed:  ',maxval(abs(u1(1:lx1))/1e3)


  !BOUNDARY CONDITIONS
  sinwt=exp(-(t-tdel)**2/tfwhm**2)*sin(omega0*t)
  sinwt2=exp(-(t-tdel)**2/(tfwhm/2)**2)*sin(2*omega0*t)

  rhom(-1:0)=rhom(1)
  rhom(lx1+1:lx1+2)=rhom(lx1)
  rhou1(-1:0)=rhou1(1)
  rhou1(lx1+1:lx1+2)=rhou1(lx1)
  rhou2(-1:0)=rhou2(1)
  rhou2(lx1+1:lx1+2)=rhou2(lx1)
  rhou3(-1:0)=rhom(0)*u31*sinwt
  rhou3(lx1+1:lx1+2)=rhom(lx1+1)*u31*sinwt2
  B2(-1:0)=B2(1)
  B2(lx1+1:lx1+2)=B2(lx1)
  B3(-1:0)=B31*sinwt
  B3(lx1+1:lx1+2)=B31*sinwt2
  rhoeps(-1:0)=rhoeps(1)
  rhoeps(lx1+1:lx1+2)=rhoeps(lx1)


  !ADVECTION SUBSTEP
  v1i(:)=0.5*(rhou1(0:lx1)/rhom(0:lx1)+rhou1(1:lx1+1)/rhom(1:lx1+1))    !CELL INTERFACE SPEEDS (LEFT WALL OF ITH CELL)

  rhom=advec1D_MC(rhom,v1i,dt,dx1,dx1i)
  rhou1=advec1D_MC(rhou1,v1i,dt,dx1,dx1i)
  rhou2=advec1D_MC(rhou2,v1i,dt,dx1,dx1i)
  rhou3=advec1D_MC(rhou3,v1i,dt,dx1,dx1i)
  B2=advec1D_MC(B2,v1i,dt,dx1,dx1i)
  B3=advec1D_MC(B3,v1i,dt,dx1,dx1i)
  rhoeps=advec1D_MC(rhoeps,v1i,dt,dx1,dx1i)


  !VON NEUMANN-RICHTMYER ARTIFICIAL VISCOSITY ('CLEANS UP' SOLUTIONS
  !WITH STRONG SHOCKS) [POTTER, 1970; STONE ET AL 1992]
  Q=0.25*xi**2*(min(du1full,0.0))**2*rhom(1:lx1)


  !SOURCE TERMS FOR 1-COMP OF MOMENTUM
  p=rhoeps(1:lx1)*(gammap-1)+(B1(1:lx1)**2+B2(1:lx1)**2+B3(1:lx1)**2)/2/mu0+Q
  rhou1(1:lx1)=rhou1(1:lx1)+dt*(-1*derivative(p,dx1(1:lx1)))


  !SOURCE TERMS FOR 2-COMP OF MOMENTUM
  rhou2(1:lx1)=rhou2(1:lx1)+dt*(1/mu0*B1(1:lx1)*derivative(B2(1:lx1),dx1(1:lx1)))


  !SOURCE TERMS FOR Y-COMP OF MOMENTUM
  rhou3(1:lx1)=rhou3(1:lx1)+dt*(1/mu0*B1(1:lx1)*derivative(B3(1:lx1),dx1(1:lx1)))


  !FLOWS FROM MOMENTA
  u1=rhou1/rhom
  u2=rhou2/rhom
  u3=rhou3/rhom


  !SOURCE TERMS FOR X-COMP OF B-FIELD
  B2(1:lx1)=B2(1:lx1)+dt*derivative(B1(1:lx1)*u2(1:lx1),dx1(1:lx1))


  !SOURCE TERMS FOR Y-COMP OF B-FIELD
  grad1B1u3=derivative(B1(1:lx1)*u3(1:lx1),dx1(1:lx1))
!  grad1B1u3(1)=(B1(2)*u3(2)-B1(1)*u3(0))/(dx1(1)+dx1(2))              !USE GHOST CELL BCS TO COUPLE WAVE PERTURBATIONS ONTO MESH (check B1(1) --> B1(0))
!  grad1B1u3(lx1)=(B1(lx1)*u3(lx1+1)-B1(lx1-1)*u3(lx1-1))/(dx1(lx1+1)+dx1(lx1))      !DITTO FOR RIGHT BOUNDARY
  grad1B1u3(1)=(B1(2)*u3(2)-B1(0)*u3(0))/(dx1(2)+dx1(1))
  grad1B1u3(lx1)=(B1(lx1+1)*u3(lx1+1)-B1(lx1-1)*u3(lx1-1))/(dx1(lx1+1)+dx1(lx1))
  B3(1:lx1)=B3(1:lx1)+dt*grad1B1u3


  !SOURCE TERMS FOR INTERNAL ENERGY (RK2 STEPPING, I THINK THIS MIGHT ACTUALLY BE DOABLE IN ONE STEP...)
  p=rhoeps(1:lx1)*(gammap-1)+Q
  grad1u1=derivative(u1(1:lx1),dx1(1:lx1))
  sourceterm = -p*grad1u1
  rhoepshalf = dt/2*sourceterm + rhoeps(1:lx1)
  p=rhoepshalf*(gammap-1)+Q
  sourceterm = -p*grad1u1
  rhoeps(1:lx1)=rhoeps(1:lx1)+dt*sourceterm


  !OUTPUT
  if (t>toutnext) then
    write(u,*) t
    call writearray(u,rhom/amu)
    call writearray(u,u1)
    call writearray(u,u2)
    call writearray(u,u3)
    call writearray(u,B1)
    call writearray(u,B2)
    call writearray(u,B3)
    call writearray(u,rhoeps*(gammap-1))
    toutnext=toutnext+dtout
  end if
end do

close(u)


contains

  subroutine writearray(fileunit,array)
    integer, intent(in) :: fileunit
    real(wp), dimension(:), intent(in) :: array

    integer :: k

    do k=1,size(array)
      write(fileunit,*) array(k)
    end do
  end subroutine writearray


  function derivative(f,dx)
    real(wp), dimension(:), intent(in) :: dx  !presumed backward diffs.
    real(wp), dimension(:), intent(in) :: f

    integer :: lx
    real(wp), dimension(1:size(f)) :: derivative

    lx=size(f)
    derivative(1)=(f(2)-f(1))/dx(2)
    derivative(2:lx-1)=(f(3:lx)-f(1:lx-2))/(dx(3:lx)+dx(2:lx-1))
    derivative(lx)=(f(lx)-f(lx-1))/dx(lx)

  end function derivative

end program MHD1D_SAW

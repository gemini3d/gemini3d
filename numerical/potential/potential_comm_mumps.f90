module potential_comm

!------------------------------------------------------------
!-------THIS MODULE IS MEANT TO WORK WITH THE MUMPS 2D
!-------INTEGRATED SOLVER IF THE GRID IS 3D, OR A FIELD-RESOLVED
!-------SOLVER IF THE GRID IS 2D (MUMPS CAN'T HANDLE 3D VERY WELL).
!------------------------------------------------------------

!NOTE THAT ONLY THE CURVILINEAR FUNCTION ARE UP TO DATE.

use mpi, only: mpi_integer, mpi_comm_world, mpi_status_ignore

use phys_consts, only : wp, pi, lsp
use grid, only : curvmesh,flagswap, gridflag
use collisions, only: conductivities, capacitance
use calculus, only: div3d, integral3d1, grad3d1, grad3d2, grad3d3, integral3d1_curv_alt
use potentialBCs_mumps, only: potentialbcs2D, potentialbcs2D_fileinput
use potential_mumps, only : elliptic3D_curv, &
                            elliptic2D_pol_conv_curv, &
                            elliptic2D_pol_conv_curv_periodic2, &
                            elliptic2D_nonint_curv, &
                            elliptic_workers

use mpimod

implicit none

private



!OVERLOAD THE SOLVERS FOR DEALING WITH THE CURVILINEAR MESHES
!NOTE WORKER SUBROUTINE DOES NOT NEED TO BE CHANGED/OVERLOADED
interface electrodynamics
  module procedure electrodynamics_curv
end interface electrodynamics

interface potential_root_mpi
  module procedure potential_root_mpi_curv
end interface potential_root_mpi


public :: electrodynamics


contains


subroutine electrodynamics_curv(it,t,dt,nn,vn2,vn3,Tn,sourcemlat,ns,Ts,vs1,B1,vs2,vs3,x, &
                         potsolve,flagcap, &
                         E1,E2,E3,J1,J2,J3,Phiall, &
                         flagE0file,dtE0,E0dir,ymd,UTsec)

!------------------------------------------------------------
!-------THIS IS A WRAPPER FUNCTION FOR THE ELECTRODYANMICS
!-------PART OF THE MODEL.  BOTH THE ROOT AND WORKER PROCESSES
!-------CALL THIS SAME SUBROUTINE, WHEN THEN BRANCHES INTO 
!-------DIFFERENT TASKS FOR EACH AFTER ALL COMPUTE CONDUCTIVITIES
!-------AND INERTIAL CAPACITANCE.
!-------
!-------NOTE THAT THE ALLOCATION STATUS
!-------OF THE ALL VARIABLES FOR THE WORKERS WILL BE UN-
!-------ALLOCATED.  SO THIS CODE IS ONLY COMPLIANT WITH 
!-------FORTRAN >= 2003 STANDARDS.
!------------------------------------------------------------

integer, intent(in) :: it
real(wp), intent(in) :: t,dt

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: vn2,vn3,Tn
real(wp), intent(in) :: sourcemlat
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  vs2,vs3

type(curvmesh), intent(in) :: x

integer, intent(in) :: potsolve
integer, intent(in) :: flagcap

real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall     !inout since it may not be allocated or deallocated in this procedure

integer, intent(in) :: flagE0file
real(wp), intent(in) :: dtE0
character(*), intent(in) :: E0dir
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec

real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: sig0,sigP,sigH
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,1:size(ns,4)) :: muP,muH,muPvn,muHvn
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: incap
real(wp) :: tstart,tfin 

integer :: lx1,lx2,lx3,isp
integer :: ix1,ix2,ix3,iinull
real(wp) :: minh1,maxh1,minh2,maxh2,minh3,maxh3


!SIZES
lx1=size(ns,1)-4
lx2=size(ns,2)-4
lx3=size(ns,3)-4


!POTENTIAL SOLUTION (IF REQUIRED)
call cpu_time(tstart)
call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,muPvn,muHvn)
if (flagcap/=0) then
  call capacitance(ns,B1,flagcap,incap)
  if (it==1) then     !check that we don't have an unsupported grid type for doing capacitance
    minh1=minval(x%h1)
    maxh1=maxval(x%h1)
    minh2=minval(x%h2)
    maxh2=maxval(x%h2)
    minh3=minval(x%h3)    !is it okay for a worker to bail on a process???
    maxh3=maxval(x%h3)
    if (minh1<0.99d0 .or. maxh1>1.01d0 .or. minh2<0.99d0 .or. maxh2>1.01d0 .or. minh3<0.99d0 .or. maxh3>1.01d0) then
      error stop 'Capacitance is being calculated for possibly unsupported grid type. Please check input file settings.'
    end if
  end if
else
  incap=0d0
end if
call cpu_time(tfin)
if (myid==0) then
  if (flagcap/=0) then
    print *, 'Conductivities and capacitance for time step:  ',t,' took ',tfin-tstart,' seconds...'
  else
    print *, 'Conductivities for time step:  ',t,' took ',tfin-tstart,' seconds...'
  end if
end if

if (potsolve == 1 .or. potsolve == 3) then    !electrostatic solve or electrostatic alt. solve
  call cpu_time(tstart)

  if (myid/=0) then    !role-specific communication pattern (all-to-root-to-all), workers initiate with sends
     call potential_workers_mpi(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,vn2,vn3,sourcemlat,B1,x, &
                            potsolve, &
                            E1,E2,E3,J1,J2,J3)
  else
    call potential_root_mpi(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,vn2,vn3,sourcemlat,B1,x, &
                              potsolve, &
                              E1,E2,E3,J1,J2,J3, &
                              Phiall,flagE0file,dtE0,E0dir,ymd,UTsec)
  end if

  !DRIFTS - NEED TO INCLUDE ELECTRIC, WIND-DRIVEN, AND GRAVITATIONAL???
  if (lx2/=1) then    !full 3D solve, go with the regular formulas
    do isp=1,lsp
      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+muPvn(:,:,:,isp)*vn2-muHvn(:,:,:,isp)*vn3
      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+muHvn(:,:,:,isp)*vn2+muPvn(:,:,:,isp)*vn3
    end do
  else                !flip signs on the cross products in 2D.  Note that due to dimension shuffling E2,3 mapping is already handled
    do isp=1,lsp
      vs2(1:lx1,1:lx2,1:lx3,isp)=-muP(:,:,:,isp)*E2+muH(:,:,:,isp)*E3-muPvn(:,:,:,isp)*vn2+muHvn(:,:,:,isp)*vn3
      vs3(1:lx1,1:lx2,1:lx3,isp)=-muH(:,:,:,isp)*E2-muP(:,:,:,isp)*E3-muHvn(:,:,:,isp)*vn2-muPvn(:,:,:,isp)*vn3
    end do
  end if

!    do isp=1,lsp    !To leading order the ion drifts do not include the polarization parts, otherwise it may mess up polarization convective term in the electrodynamics solver...
!      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE2Dt
!      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE3Dt
!    end do

  call cpu_time(tfin)

  if (myid==0) then
    print *, 'Potential solution for time step:  ',t,' took ',tfin-tstart,' seconds...'
    print *, 'Min and max root drift values:  ',minval(vs2),maxval(vs2), minval(vs3),maxval(vs3)
  end if

else if (potsolve == 2) then  !inductive form of model, could this be subcycled to speed things up?
  !Do nothing for now...
else   !null solve; just force everything to zero
  E1=0d0; E2=0d0; E3=0d0; J1=0d0; J2=0d0; J3=0d0;
  vs2=0d0; vs3=0d0;
end if

end subroutine electrodynamics_curv


subroutine potential_root_mpi_curv(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,vn2,vn3,sourcemlat,B1,x, &
                            potsolve, &
                            E1,E2,E3,J1,J2,J3, &
                            Phiall,flagE0file,dtE0,E0dir,ymd,UTsec)

!------------------------------------------------------------
!-------ROOT MPI COMM./SOLVE ROUTINE FOR POTENTIAL.  THIS VERSION
!-------INCLUDES THE POLARIZATION CURRENT TIME DERIVATIVE PART
!-------AND CONVECTIVE PARTS IN MATRIX SOLUTION.
!-------STATE VARIABLES VS2,3 INCLUDE GHOST CELLS.  FOR NOW THE
!-------POLARIZATION TERMS ARE PASSED BACK TO MAIN FN, EVEN THOUGH
!-------THEY ARE NOT USED (THEY MAY BE IN THE FUTURE)
!------------------------------------------------------------

integer, intent(in) :: it
real(wp), intent(in) :: t,dt
real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH
real(wp), dimension(:,:,:), intent(in) ::  incap
real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3    
real(wp), dimension(:,:,:), intent(in) ::  vn2,vn3
real(wp), intent(in) :: sourcemlat
real(wp), dimension(-1:,-1:,-1:), intent(in) ::  B1

type(curvmesh), intent(in) :: x

integer, intent(in) :: potsolve

real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
real(wp), dimension(:,:,:), intent(inout) :: Phiall   !not good form, but I'm lazy...  Forgot what I meant by this...

integer, intent(in) :: flagE0file
real(wp), intent(in) :: dtE0
character(*), intent(in) :: E0dir
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: v2,v3 

real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: srctermall
real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)), target :: Vminx1,Vmaxx1     !allow pointer aliases for these vars.
real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)) :: Vminx1buf,Vmaxx1buf
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,3)) :: Vminx2,Vmaxx2
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2)) :: Vminx3,Vmaxx3
integer :: flagdirich

real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)) :: v2slaball,v3slaball   !stores drift velocs. for pol. current

!   real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: paramtrim    !to hold trimmed magnetic field

real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: E01all,E02all,E03all    !background fields
!   real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: divJperpall2,divJperpall3,divJperpall   !more work arrays

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: integrand,sigintegral    !general work array for doing integrals

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: grad2E,grad3E    !more work arrays for pol. curr.
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: DE2Dt,DE3Dt   !pol. drift
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: J1pol,J2pol,J3pol

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E01,E02,E03   !distributed background fields
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: srcterm,divJperp
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E1prev,E2prev,E3prev
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: Phi

real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: SigPint2,SigPint3,SigHint,incapint,srctermint
real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)) :: SigPint2all,SigPint3all,SigHintall,incapintall,srctermintall

real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)) :: Phislab,Phislab0

real(wp), dimension(0:size(E1,1)+1,0:size(E1,2)+1,0:size(E1,3)+1) :: divtmp    !one extra grid point on either end to facilitate derivatives
real(wp), dimension(-1:size(E1,1)+2,-1:size(E1,2)+2,-1:size(E1,3)+2) :: J1halo,J2halo,J3halo    !haloing assumes existence of two ghost cells

real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: sig0scaledall,sigPscaledall,sigHscaledall
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: sig0scaled,sigPscaled,sigHscaled

logical :: perflag    !MUMPS stuff

real(wp), dimension(1:size(Phiall,3)) :: Vminx2slice,Vmaxx2slice
real(wp), dimension(1:size(Phiall,2)) :: Vminx3slice,Vmaxx3slice
real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: Vminx1slab,Vmaxx1slab
real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: v2slab,v3slab

integer :: iid
integer :: ix1,ix2,ix3,lx1,lx2,lx3,lx3all
integer :: idleft,idright,iddown,idup

real(wp) :: tstart,tfin


!SIZES - PERHAPS SHOULD BE TAKEN FROM GRID MODULE INSTEAD OF RECOMPUTED?
lx1=size(sig0,1)
lx2=size(sig0,2)
lx3=size(sig0,3)
lx3all=size(Phiall,3)


!USE PREVIOUS MUMPS PERMUTATION (OLD CODE? BUT MIGHT BE WORTH REINSTATING?)
!    perflag=.false.
perflag=.true.


!R-------
!POPULATE BACKGROUND AND BOUNDARY CONDITION ARRAYS - IDEALLY ROOT ONLY SINCE IT INVOLVES FILE INPUT, ALTHOUGH THE INTERPOLATION MAY BE SLOW...
call cpu_time(tstart)
if (flagE0file==1) then
  call potentialBCs2D_fileinput(dt,dtE0,t,ymd,UTsec,E0dir,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                      E01all,E02all,E03all,flagdirich) 
else
  call potentialBCs2D(t,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                      E01all,E02all,E03all,flagdirich)     !user needs to manually swap x2 and x3 in this function.
end if
call cpu_time(tfin)
print *, 'Root has computed BCs in time:  ',tfin-tstart
!R-------


!R--------
do iid=1,lid-1    !communicate intent for solve to workers so they know whether or not to call mumps fn.
  call mpi_send(flagdirich,1,MPI_INTEGER,iid,tagflagdirich,MPI_COMM_WORLD,ierr)
end do
print *, 'Root has communicated type of solve to workers:  ',flagdirich
!R--------


!Need to broadcast background fields from root
!Need to also broadcast x1 boundary conditions for source term calculations.  
call bcast_send(E01all,tagE01,E01)
call bcast_send(E02all,tagE02,E02)
call bcast_send(E03all,tagE03,E03)

!These are pointer targets so don't assume contiguous in memory - pack them into a buffer to be safe
Vmaxx1buf=Vmaxx1; Vminx1buf=Vminx1;
call bcast_send(Vminx1buf,tagVminx1,Vminx1slab)
call bcast_send(Vmaxx1buf,tagVmaxx1,Vmaxx1slab)


!-------
!CONDUCTION CURRENT BACKGROUND SOURCE TERMS FOR POTENTIAL EQUATION. MUST COME AFTER CALL TO BC CODE.
J1=0d0    !so this div is only perp components
if (flagswap==1) then
  J2=sigP*E02+sigH*E03    !BG x2 current
  J3=-1*sigH*E02+sigP*E03    !BG x3 current
else
  J2=sigP*E02-sigH*E03    !BG x2 current
  J3=sigH*E02+sigP*E03    !BG x3 current
end if
!    srcterm=div3D(J1,J2,J3,x,1,lx1,1,lx2,1,lx3)     !first part of source term so do a straight assignment here
J1halo(1:lx1,1:lx2,1:lx3)=J1
J2halo(1:lx1,1:lx2,1:lx3)=J2
J3halo(1:lx1,1:lx2,1:lx3)=J3

J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)

J2halo(0,1:lx2,1:lx3)=J2halo(1,1:lx2,1:lx3)
J2halo(lx1+1,1:lx2,1:lx3)=J2halo(lx1,1:lx2,1:lx3)

J3halo(0,1:lx2,1:lx3)=J3halo(1,1:lx2,1:lx3)
J3halo(lx1+1,1:lx2,1:lx3)=J3halo(lx1,1:lx2,1:lx3)

call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 
call halo(J2halo,1,tagJ2,x%flagper)
call halo(J3halo,1,tagJ3,x%flagper)

!ZZZ - NEED TO SET GLOBAL BOUNDARY FOR APPROPRIATE WORKERS USING ZOH EXTRAP. (unless periodic)
idleft=myid3-1
idright=myid3+1
iddown=myid2-1
idup=myid2+1
if (iddown==-1) then
  J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  J2halo(1:lx1,0,1:lx3)=J2halo(1:lx1,1,1:lx3)
  J3halo(1:lx1,0,1:lx3)=J3halo(1:lx1,1,1:lx3)
end if
if (idup==lid2) then
  J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
  J2halo(1:lx1,lx2+1,1:lx3)=J2halo(1:lx1,lx2,1:lx3)
  J3halo(1:lx1,lx2+1,1:lx3)=J3halo(1:lx1,lx2,1:lx3)
end if
if (.not. x%flagper) then
  if (idleft==-1) then
    J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
    J2halo(1:lx1,1:lx2,0)=J2halo(1:lx1,1:lx2,1)
    J3halo(1:lx1,1:lx2,0)=J3halo(1:lx1,1:lx2,1)
  end if
  if (idright==lid3) then
    J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
    J2halo(1:lx1,1:lx2,lx3+1)=J2halo(1:lx1,1:lx2,lx3)
    J3halo(1:lx1,1:lx2,lx3+1)=J3halo(1:lx1,1:lx2,lx3)
  end if
end if

divtmp=div3D(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),J2halo(0:lx1+1,0:lx2+1,0:lx3+1), &
             J3halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
srcterm=divtmp(1:lx1,1:lx2,1:lx3)
print *, 'Root has computed background field source terms...',minval(srcterm), maxval(srcterm)
!-------


!-------
!NEUTRAL WIND SOURCE TERMS FOR POTENTIAL EQUATION, SIMILAR TO ABOVE BLOCK OF CODE
J1=0d0    !so this div is only perp components
if (flagswap==1) then
  J2=-1*sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)     !wind x2 current, note that all workers already have a copy of this.  
  J3=sigH*vn3*B1(1:lx1,1:lx2,1:lx3)+sigP*vn2*B1(1:lx1,1:lx2,1:lx3)     !wind x3 current
else
  J2=sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)     !wind x2 current
  J3=sigH*vn3*B1(1:lx1,1:lx2,1:lx3)-sigP*vn2*B1(1:lx1,1:lx2,1:lx3)     !wind x3 current
end if
!    divJperp=div3D(J1,J2,J3,x,1,lx1,1,lx2,1,lx3)
J1halo(1:lx1,1:lx2,1:lx3)=J1
J2halo(1:lx1,1:lx2,1:lx3)=J2
J3halo(1:lx1,1:lx2,1:lx3)=J3

J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)

J2halo(0,1:lx2,1:lx3)=J2halo(1,1:lx2,1:lx3)
J2halo(lx1+1,1:lx2,1:lx3)=J2halo(lx1,1:lx2,1:lx3)

J3halo(0,1:lx2,1:lx3)=J3halo(1,1:lx2,1:lx3)
J3halo(lx1+1,1:lx2,1:lx3)=J3halo(lx1,1:lx2,1:lx3)

call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 
call halo(J2halo,1,tagJ2,x%flagper)
call halo(J3halo,1,tagJ3,x%flagper)

!CHECK AND SET GLOBAL BOUNDARY AFTER HALO OPERATIONS
if (iddown==-1) then
  J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  J2halo(1:lx1,0,1:lx3)=J2halo(1:lx1,1,1:lx3)
  J3halo(1:lx1,0,1:lx3)=J3halo(1:lx1,1,1:lx3)
end if
if (idup==lid2) then
  J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
  J2halo(1:lx1,lx2+1,1:lx3)=J2halo(1:lx1,lx2,1:lx3)
  J3halo(1:lx1,lx2+1,1:lx3)=J3halo(1:lx1,lx2,1:lx3)
end if
if (.not. x%flagper) then
  if (idleft==-1) then
    J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
    J2halo(1:lx1,1:lx2,0)=J2halo(1:lx1,1:lx2,1)
    J3halo(1:lx1,1:lx2,0)=J3halo(1:lx1,1:lx2,1)
  end if
  if (idright==lid3) then
    J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
    J2halo(1:lx1,1:lx2,lx3+1)=J2halo(1:lx1,1:lx2,lx3)
    J3halo(1:lx1,1:lx2,lx3+1)=J3halo(1:lx1,1:lx2,lx3)
  end if
end if

divtmp=div3D(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),J2halo(0:lx1+1,0:lx2+1,0:lx3+1), &
             J3halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
srcterm=srcterm+divtmp(1:lx1,1:lx2,1:lx3)
print *, 'Root has computed wind source terms...',minval(srcterm),  maxval(srcterm)
!-------


!!!!!!!!
!-----AT THIS POINT WE MUST DECIDE WHETHER TO DO AN INTEGRATED SOLVE OR A 2D FIELD-RESOLVED SOLVED
!-----DECIDE BASED ON THE SIZE OF THE X2 DIMENSION
if (lx2/=1) then    !either field-resolved 3D or integrated 2D solve for 3D domain
  if (potsolve == 1) then    !2D, field-integrated solve
    print *, 'Beginning field-integrated solve...'


    !> INTEGRATE CONDUCTANCES AND CAPACITANCES FOR SOLVER COEFFICIENTS 
    integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)
    sigintegral=integral3D1(integrand,x,1,lx1)    !no haloing required for a field-line integration
    SigPint2=sigintegral(lx1,:,:)

    integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)
    sigintegral=integral3D1(integrand,x,1,lx1)
    SigPint3=sigintegral(lx1,:,:)

    integrand=x%h1(1:lx1,1:lx2,1:lx3)*sigH
    sigintegral=integral3D1(integrand,x,1,lx1)
    SigHint=sigintegral(lx1,:,:)
    
    sigintegral=integral3D1(incap,x,1,lx1)
    incapint=sigintegral(lx1,:,:)
    !-------


    !> PRODUCE A FIELD-INTEGRATED SOURCE TERM
    if (flagdirich /= 1) then   !Neumann conditions; incorporate a source term and execute the solve
      print *, 'Using FAC boundary condition...'
      !-------
      integrand=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)*srcterm
      sigintegral=integral3D1(integrand,x,1,lx1)
      srctermint=sigintegral(lx1,:,:)
      srctermint=srctermint+x%h2(lx1,1:lx2,1:lx3)*x%h3(lx1,1:lx2,1:lx3)*Vmaxx1slab- &
                                 x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab     !workers don't have access to boundary conditions, unless root sends
      !-------


      v2slab=vs2(lx1,1:lx2,1:lx3,1); v3slab=vs3(lx1,1:lx2,1:lx3,1);    !need to pick out the ExB drift here (i.e. the drifts from highest altitudes); but this is only valid for Cartesian, so it's okay for the foreseeable future


      !RADD--- ROOT NEEDS TO PICK UP *INTEGRATED* SOURCE TERMS AND COEFFICIENTS FROM WORKERS
      call gather_recv(srctermint,tagsrc,srctermintall)
      call gather_recv(incapint,tagincapint,incapintall)
      call gather_recv(SigPint2,tagSigPint2,SigPint2all)
      call gather_recv(SigPint3,tagSigPint3,SigPint3all)
      call gather_recv(SigHint,tagSigHint,SigHintall)
      call gather_recv(v2slab,tagv2electro,v2slaball)
      call gather_recv(v3slab,tagv3electro,v3slaball)


      !R------
      !EXECUTE FIELD-INTEGRATED SOLVE
      Vminx2slice=Vminx2(lx1,:)    !slice the boundaries into expected shape
      Vmaxx2slice=Vmaxx2(lx1,:)
      Vminx3slice=Vminx3(lx1,:)
      Vmaxx3slice=Vmaxx3(lx1,:)
      Phislab0=Phiall(lx1,:,:)    !root already possess the fullgrid potential from prior solves...
      print *, 'Root is calling MUMPS...'
      !R-------

      !R------ EXECUTE THE MUMPS SOLVE FOR FIELD-INT
      call cpu_time(tstart)
      if (.not. x%flagper) then     !nonperiodic mesh
        print *, '!!!User selected aperiodic solve...'
        Phislab=elliptic2D_pol_conv_curv(srctermintall,SigPint2all,SigPint3all,SigHintall,incapintall,v2slaball,v3slaball, &
                                 Vminx2slice,Vmaxx2slice,Vminx3slice,Vmaxx3slice, &
                                 dt,x,Phislab0,perflag,it)    !note tha this solver is only valid for cartesian meshes, unless the inertial capacitance is set to zero
      else
        print *, '!!!User selected periodic solve...'
        Phislab=elliptic2D_pol_conv_curv_periodic2(srctermintall,SigPint2all,SigHintall,incapintall,v2slaball,v3slaball, &    !note that either sigPint2 or 3 will work since this must be cartesian...
                                   Vminx2slice,Vmaxx2slice,Vminx3slice,Vmaxx3slice, &
                                   dt,x,Phislab0,perflag,it)
      end if
      call cpu_time(tfin)
      print *, 'Root received results from MUMPS which took time:  ',tfin-tstart
      !R-------

    else     !Dirichlet conditions - since this is field integrated we just copy BCs specified by user to other locations along field line
      !R------
      Phislab=Vmaxx1   !potential is whatever user specifies, since we assume equipotential field lines, it doesn't really matter whether we use Vmaxx1 or Vminx1.  Note however, tha thte boundary conditions subroutines should explicitly set these to be equal with Dirichlet conditions, for consistency.
      print *, 'Dirichlet conditions selected with field-integrated solve. Copying BCs along x1-direction...'
      !R------
    end if


    !R------  AFTER ANY TYPE OF FIELD-INT SOLVE COPY THE BCS ACROSS X1 DIMENSION
    do ix1=1,lx1
      Phiall(ix1,:,:)=Phislab(:,:) !copy the potential across the ix1 direction; past this point there is no real difference with 3D, note that this is still valid in curvilinear form
    end do
    !R------

  else    !resolved 3D solve
    !ZZZ - conductivities need to be properly scaled here...  So does the source term...  Maybe leave as broken for now since I don't really plan to use this code


    !RADD--- ROOT NEEDS TO PICK UP FIELD-RESOLVED SOURCE TERM AND COEFFICIENTS FROM WORKERS
    call gather_recv(sigPscaled,tagsigP,sigPscaledall)
    call gather_recv(sigHscaled,tagsigH,sigHscaledall)
    call gather_recv(sig0scaled,tagsig0,sig0scaledall)
    call gather_recv(srcterm,tagsrc,srctermall)

    !R------
    print *, '!Beginning field-resolved 3D solve (could take a very long time)...'
    Phiall=elliptic3D_curv(srctermall,sig0scaledall,sigPscaledall,sigHscaledall,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                      x,flagdirich,perflag,it)
    !R------
  end if
else   !lx1=1 so do a field-resolved 2D solve over x1,x3
  print *, 'Beginning field-resolved 2D solve...  Type;  ',flagdirich


  !-------
  !PRODUCE SCALED CONDUCTIVITIES TO PASS TO SOLVER, ALSO SCALED SOURCE TERM
  sig0scaled=x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)*sig0
  if (flagswap==1) then
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)*sigP    !remember to swap 2-->3
  else
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
  end if
  srcterm=srcterm*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
!      srcterm=-1d0*srcterm     !in a 2D solve negate this due to it being a cross produce and the fact that we've permuted the 2 and 3 dimensions.  ZZZ - NOT JUST THIS WORKS WITH BACKGROUND FIELDS???
  !-------

  !RADD--- NEED TO GET THE RESOLVED SOURCE TERMS AND COEFFICIENTS FROM WORKERS
  call gather_recv(sigPscaled,tagsigP,sigPscaledall)
  call gather_recv(sig0scaled,tagsig0,sig0scaledall)
  call gather_recv(srcterm,tagsrc,srctermall)


  !EXECUTE THE SOLVE WITH MUMPS AND SCALED TERMS
  !NOTE THE LACK OF A SPECIAL CASE HERE TO CHANGE THE POTENTIAL PROBLEM - ONLY THE HALL TERM CHANGES (SINCE RELATED TO EXB) BUT THAT DOESN'T APPEAR IN THIS EQN!
  Phiall=elliptic2D_nonint_curv(srctermall,sig0scaledall,sigPscaledall,Vminx1,Vmaxx1,Vminx3,Vmaxx3, &
                    x,flagdirich,perflag,it)
end if
print *, 'MUMPS time:  ',tfin-tstart
!!!!!!!!!


!RADD--- ROOT NEEDS TO PUSH THE POTENTIAL BACK TO ALL WORKERS FOR FURTHER PROCESSING (BELOW)
call bcast_send(Phiall,tagPhi,Phi)


!-------
!STORE PREVIOUS TIME TOTAL FIELDS BEFORE UPDATING THE ELECTRIC FIELDS WITH NEW POTENTIAL (OLD FIELDS USED TO CALCULATE POLARIZATION CURRENT)
E1prev=E1
E2prev=E2
E3prev=E3
!-------


!-------
!CALCULATE PERP FIELDS FROM POTENTIAL
!      E20all=grad3D2(-1d0*Phi0all,dx2(1:lx2))    !causes major memory leak. maybe from arithmetic statement argument? Left here as a 'lesson learned' (or is it a gfortran bug...)
!      E30all=grad3D3(-1d0*Phi0all,dx3all(1:lx3all))
Phi=-1d0*Phi
E2=grad3D2(Phi,x,1,lx1,1,lx2,1,lx3)    !no haloing required
!    E3=grad3D3(Phi,x,1,lx1,1,lx2,1,lx3)    !needs to be haloed

J1halo(1:lx1,1:lx2,1:lx3)=Phi
J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)
call halo(J1halo,1,tagJ1,x%flagper)

if (iddown==-1) then
  J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
end if
if (idup==lid2) then
  J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
end if
if (.not. x%flagper) then
  if (idleft==-1) then
    J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
  end if
  if (idright==lid3) then
    J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
  end if
end if

divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
E3=divtmp(1:lx1,1:lx2,1:lx3)
Phi=-1d0*Phi   !put things back for later use
!--------


!R-------
!JUST TO JUDGE THE IMPACT OF MI COUPLING
print *, 'Max integrated inertial capacitance:  ',maxval(incapint)
print *, 'Max integrated Pedersen conductance (includes metric factors):  ',maxval(SigPint2)
print *, 'Max integrated Hall conductance (includes metric factors):  ',minval(SigHint), maxval(SigHint)
print *, 'Max E2,3 BG and response values are:  ',maxval(E02), maxval(E03),maxval(E2),maxval(E3)
print *, 'Min E2,3 BG and response values are:  ',minval(E02), minval(E03),minval(E2),minval(E3)
!R-------


!--------
!ADD IN BACKGROUND FIELDS BEFORE DISTRIBUTING TO WORKERS
E2=E2+E02
E3=E3+E03
!--------


!--------
!COMPUTE TIME DERIVATIVE NEEDED FOR POLARIZATION CURRENT.  ONLY DO THIS IF WE HAVE SPECIFIC NONZERO INERTIAL CAPACITANCE
if (maxval(incap) > 1d-6) then
  grad2E=grad3D2(E2,x,1,lx1,1,lx2,1,lx3)
!      grad3E=grad3D3(E2,x,1,lx1,1,lx2,1,lx3)

  J1halo(1:lx1,1:lx2,1:lx3)=E2
  J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
  J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)

  call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 

  if (iddown==-1) then
    J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  end if
  if (idup==lid2) then
    J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
  end if
  if (.not. x%flagper) then
    if (idleft==-1) then
      J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
    end if
    if (idright==lid3) then
      J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
    end if
  end if
  
  divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad3E=divtmp(1:lx1,1:lx2,1:lx3)

  DE2Dt=(E2-E2prev)/dt+v2*grad2E+v3*grad3E
  grad2E=grad3D2(E3,x,1,lx1,1,lx2,1,lx3)
!      grad3E=grad3D3(E3,x,1,lx1,1,lx2,1,lx3)

  J1halo(1:lx1,1:lx2,1:lx3)=E3
  J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
  J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)
  call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 

  if (iddown==-1) then
    J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  end if
  if (idup==lid2) then
    J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
  end if
  if (.not. x%flagper) then
    if (idleft==-1) then
      J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
    end if
    if (idright==lid3) then
      J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
    end if
  end if

  divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad3E=divtmp(1:lx1,1:lx2,1:lx3)

  DE3Dt=(E3-E3prev)/dt+v2*grad2E+v3*grad3E
  J1pol=0d0
  J2pol=incap*DE2Dt
  J3pol=incap*DE3Dt
else       !pure electrostatic solve was done
  DE2Dt=0d0
  DE3Dt=0d0
  J1pol=0d0 
  J2pol=0d0
  J3pol=0d0      
end if
!--------

!-------
if (flagswap==1) then
  J2=sigP*E2+sigH*E3    !BG field already added to E above
  J3=-1*sigH*E2+sigP*E3
else
  J2=sigP*E2-sigH*E3    !BG field already added to E above
  J3=sigH*E2+sigP*E3
end if

!WHAT I THINK THE NEUTRAL WIND CURRENTS SHOULD BE IN 2D
if (flagswap==1) then
  J2=J2-sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
  J3=J3+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)+sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
else
  J2=J2+sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
  J3=J3+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)-sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
end if
!-------


!!!!!!!!
!NOW DEAL WITH THE PARALLEL FIELDS AND ALL CURRENTS
if (lx2/=1 .and. potsolve ==1) then    !we did a field-integrated solve above

  !-------
  !NOTE THAT A DIRECT E1ALL CALCULATION WILL GIVE ZERO, SO USE INDIRECT METHOD, AS FOLLOWS
  J1=0d0    !a placeholder so that only the perp divergence is calculated - will get overwritten later.
!      divJperp=div3D(J1,J2,J3,x,1,lx1,1,lx2,1,lx3)
  J1halo(1:lx1,1:lx2,1:lx3)=J1
  J2halo(1:lx1,1:lx2,1:lx3)=J2
  J3halo(1:lx1,1:lx2,1:lx3)=J3

  J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
  J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)
  J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)

  J2halo(0,1:lx2,1:lx3)=J2halo(1,1:lx2,1:lx3)
  J2halo(lx1+1,1:lx2,1:lx3)=J2halo(lx1,1:lx2,1:lx3)
  J2halo(1:lx1,0,1:lx3)=J2halo(1:lx1,1,1:lx3)
  J2halo(1:lx1,lx2+1,1:lx3)=J2halo(1:lx1,lx2,1:lx3)

  J3halo(0,1:lx2,1:lx3)=J3halo(1,1:lx2,1:lx3)
  J3halo(lx1+1,1:lx2,1:lx3)=J3halo(lx1,1:lx2,1:lx3)
  J3halo(1:lx1,0,1:lx3)=J3halo(1:lx1,1,1:lx3)
  J3halo(1:lx1,lx2+1,1:lx3)=J3halo(1:lx1,lx2,1:lx3)

  call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 
  call halo(J2halo,1,tagJ2,x%flagper)
  call halo(J3halo,1,tagJ3,x%flagper)

  !ZZZ - NEED TO SET GLOBAL BOUNDARIES HERE
  if (iddown==-1) then
    J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
    J2halo(1:lx1,0,1:lx3)=J2halo(1:lx1,1,1:lx3)
    J3halo(1:lx1,0,1:lx3)=J3halo(1:lx1,1,1:lx3)
  end if
  if (idup==lid2) then
    J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
    J2halo(1:lx1,lx2+1,1:lx3)=J2halo(1:lx1,lx2,1:lx3)
    J3halo(1:lx1,lx2+1,1:lx3)=J3halo(1:lx1,lx2,1:lx3)
  end if
  if (.not. x%flagper) then
    if (idleft==-1) then
      J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
      J2halo(1:lx1,1:lx2,0)=J2halo(1:lx1,1:lx2,1)
      J3halo(1:lx1,1:lx2,0)=J3halo(1:lx1,1:lx2,1)
    end if
    if (idright==lid3) then
      J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
      J2halo(1:lx1,1:lx2,lx3+1)=J2halo(1:lx1,1:lx2,lx3)
      J3halo(1:lx1,1:lx2,lx3+1)=J3halo(1:lx1,1:lx2,lx3)
    end if
  end if
  
  divtmp=div3D(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),J2halo(0:lx1+1,0:lx2+1,0:lx3+1), &
               J3halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  divJperp=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)*divtmp(1:lx1,1:lx2,1:lx3)
  if (flagdirich /= 1) then     !Neumann conditions, this is boundary location-agnostic since both bottom and top FACs are known - they have to  be loaded into VVmaxx1 and Vminx1
    print *, 'Nuemann boundaries, integrated from highest altitude down to preserve accuracy...'
    if (gridflag==0) then     !closed dipole grid, really would be best off integrating from the source hemisphere
      print *,  'Closed dipole grid; integration starting in source hemisphere (if applicable)...', &
                     minval(Vmaxx1slab), &
                     maxval(Vmaxx1slab)
      if (sourcemlat>=0d0) then    !integrate from northern hemisphere
        print *, 'Source is in northern hemisphere (or there is no source)...'
        J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,:,:))
        end do
      else
        print *, 'Source in southern hemisphere...'
        J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current starting from minx1
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,:,:))
        end do
      end if
    elseif (gridflag==1) then    !this would be an inverted grid, this max altitude corresponds to the min value of x1
      print *,  'Inverted grid; integration starting at min x1 (highest alt. or southern hemisphere)...', &
                     minval(Vminx1slab), & 
                     maxval(Vminx1slab)
      J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
      do ix1=1,lx1
        J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                         (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,:,:))
      end do
    else        !minx1 is at teh bottom of the grid to integrate from max x1
      print *,  'Non-inverted grid; integration starting at max x1...', minval(Vmaxx1slab), maxval(Vmaxx1slab)
      J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
      do ix1=1,lx1
        J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                         (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,:,:))
      end do
    end if
  else      !Dirichlet conditions - we need to integrate from the ***lowest altitude*** (where FAC is known to be zero, note this is not necessarilty the logical bottom of the grid), upwards (to where it isn't)
    if (gridflag/=2) then    !inverted grid (logical top is the lowest altitude)
      print *, 'Inverted grid detected - integrating logical top downward to compute FAC...'
      J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current
      do ix1=1,lx1
        J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                         (J1(ix1,:,:))    !FAC AT TOP ASSUMED TO BE ZERO
      end do
    else      !non-inverted grid (logical bottom is the lowest altitude - so integrate normy)
      print *, 'Non-inverted grid detected - integrating logical bottom to top to compute FAC...'
      J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
      do ix1=1,lx1
        J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                         (-1d0*J1(ix1,:,:))    !FAC AT THE BOTTOM ASSUMED TO BE ZERO
      end do
    end if
  end if
  E1=J1/sig0
  !-------

else   !we resolved the field line (either 2D solve or full 3D) so just differentiate normally

  !-------
  Phi=-1d0*Phi
  E1=grad3D1(Phi,x,1,lx1,1,lx2,1,lx3)    !no haloing required since x1-derivative
  Phi=-1d0*Phi
  J1=sig0*E1
  !-------

end if
!!!!!!!!!


!R-------
print *, 'Max topside FAC (abs. val.) computed to be:  ',maxval(abs(J1(1,:,:)))    !ZZZ - this rey needsz to be current at the "top"
print *, 'Max polarization J2,3 (abs. val.) computed to be:  ',maxval(abs(J2pol)), &
             maxval(abs(J3pol))
!    print *, 'Max conduction J2,3 (abs. val.) computed to be:  ',maxval(abs(J2)), &
!                 maxval(abs(J3))
print *, 'Max conduction J2,3  computed to be:  ',maxval(J2), &
             maxval(J3)
print *, 'Min conduction J2,3  computed to be:  ',minval(J2), &
             minval(J3)
print *, 'Max conduction J1 (abs. val.) computed to be:  ',maxval(abs(J1))
print *, 'flagswap:  ',flagswap
!R-------


!-------
!GRAND TOTAL FOR THE CURRENT DENSITY:  TOSS IN POLARIZATION CURRENT SO THAT OUTPUT FILES ARE CONSISTENT
J1=J1+J1pol
J2=J2+J2pol
J3=J3+J3pol
!-------

end subroutine potential_root_mpi_curv


subroutine potential_workers_mpi(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,vn2,vn3,sourcemlat,B1,x, &
                            potsolve, &
                            E1,E2,E3,J1,J2,J3)

!------------------------------------------------------------
!-------ROOT MPI COMM./SOLVE ROUTINE FOR POTENTIAL.  THIS VERSION
!-------INCLUDES THE POLARIZATION CURRENT TIME DERIVATIVE PART
!-------AND CONVECTIVE PARTS IN MATRIX SOLUTION.
!-------STATE VARIABLES VS2,3 INCLUDE GHOST CELLS.  FOR NOW THE
!-------POLARIZATION TERMS ARE PASSED BACK TO MAIN FN, EVEN THOUGH
!-------THEY ARE NOT USED (THEY MAY BE IN THE FUTURE)
!------------------------------------------------------------

integer, intent(in) :: it
real(wp), intent(in) :: t,dt
real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH
real(wp), dimension(:,:,:), intent(in) ::  incap
real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3    
real(wp), dimension(:,:,:), intent(in) ::  vn2,vn3
real(wp), intent(in) :: sourcemlat
real(wp), dimension(-1:,-1:,-1:), intent(in) ::  B1

type(curvmesh), intent(in) :: x

integer, intent(in) :: potsolve

real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3

integer :: flagdirich

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: paramtrim    !to hold trimmed magnetic field

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: grad2E,grad3E    !more work arrays for pol. curr.
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: DE2Dt,DE3Dt   !pol. drift
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: J1pol,J2pol,J3pol

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E01,E02,E03   !distributed background fields
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: srcterm,divJperp
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E1prev,E2prev,E3prev
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: Phi

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: integrand,sigintegral    !general work array for doing integrals
real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: SigPint2,SigPint3,SigHint,incapint,srctermint

real(wp), dimension(0:size(E1,1)+1,0:size(E1,2)+1,0:size(E1,3)+1) :: divtmp    !one extra grid point on either end to facilitate derivatives
real(wp), dimension(-1:size(E1,1)+2,-1:size(E1,2)+2,-1:size(E1,3)+2) :: J1halo,J2halo,J3halo    !haloing assumes existence of two ghost cells

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: sig0scaled,sigPscaled,sigHscaled

logical :: perflag    !MUMPS stuff

real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: Vminx1slab,Vmaxx1slab

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: v2,v3
real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: v2slab,v3slab

integer :: ix1,ix2,ix3,lx1,lx2,lx3,lx3all
integer :: idleft,idright,iddown,idup

real(wp) :: tstart,tfin


!SIZES - PERHAPS SHOULD BE TAKEN FROM GRID MODULE INSTEAD OF RECOMPUTED?
lx1=size(sig0,1)
lx2=size(sig0,2)
lx3=size(sig0,3)


!USE PREVIOUS MUMPS PERMUTATION (OLD CODE? BUT MIGHT BE WORTH REINSTATING?)
perflag=.false.


call mpi_recv(flagdirich,1,MPI_INTEGER,0,tagflagdirich,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)


!Need to broadcast background fields from root
!Need to also broadcast x1 boundary conditions for source term calculations.  
call bcast_recv(E01,tagE01)
call bcast_recv(E02,tagE02)
call bcast_recv(E03,tagE03)
call bcast_recv(Vminx1slab,tagVminx1)
call bcast_recv(Vmaxx1slab,tagVmaxx1)


!-------
!CONDUCTION CURRENT BACKGROUND SOURCE TERMS FOR POTENTIAL EQUATION. MUST COME AFTER CALL TO BC CODE.
J1=0d0    !so this div is only perp components
if (flagswap==1) then
  J2=sigP*E02+sigH*E03    !BG x2 current
  J3=-1*sigH*E02+sigP*E03    !BG x3 current
else
  J2=sigP*E02-sigH*E03    !BG x2 current
  J3=sigH*E02+sigP*E03    !BG x3 current
end if
!    srcterm=div3D(J1,J2,J3,x,1,lx1,1,lx2,1,lx3)     !first part of source term so do a straight assignment here
J1halo(1:lx1,1:lx2,1:lx3)=J1
J2halo(1:lx1,1:lx2,1:lx3)=J2
J3halo(1:lx1,1:lx2,1:lx3)=J3

J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)

J2halo(0,1:lx2,1:lx3)=J2halo(1,1:lx2,1:lx3)
J2halo(lx1+1,1:lx2,1:lx3)=J2halo(lx1,1:lx2,1:lx3)

J3halo(0,1:lx2,1:lx3)=J3halo(1,1:lx2,1:lx3)
J3halo(lx1+1,1:lx2,1:lx3)=J3halo(lx1,1:lx2,1:lx3)

call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 
call halo(J2halo,1,tagJ2,x%flagper)
call halo(J3halo,1,tagJ3,x%flagper)

!ZZZ - CHECK FOR AND SET GLOBAL BOUNDARIES
idleft=myid3-1
idright=myid3+1
iddown=myid2-1
idup=myid2+1
if (iddown==-1) then
  J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  J2halo(1:lx1,0,1:lx3)=J2halo(1:lx1,1,1:lx3)
  J3halo(1:lx1,0,1:lx3)=J3halo(1:lx1,1,1:lx3)
end if
if (idup==lid2) then
  J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
  J2halo(1:lx1,lx2+1,1:lx3)=J2halo(1:lx1,lx2,1:lx3)
  J3halo(1:lx1,lx2+1,1:lx3)=J3halo(1:lx1,lx2,1:lx3)
end if
if (.not. x%flagper) then
  if (idleft==-1) then
    J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
    J2halo(1:lx1,1:lx2,0)=J2halo(1:lx1,1:lx2,1)
    J3halo(1:lx1,1:lx2,0)=J3halo(1:lx1,1:lx2,1)
  end if
  if (idright==lid3) then
    J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
    J2halo(1:lx1,1:lx2,lx3+1)=J2halo(1:lx1,1:lx2,lx3)
    J3halo(1:lx1,1:lx2,lx3+1)=J3halo(1:lx1,1:lx2,lx3)
  end if
end if

divtmp=div3D(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),J2halo(0:lx1+1,0:lx2+1,0:lx3+1), &
             J3halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
srcterm=divtmp(1:lx1,1:lx2,1:lx3)
!-------


!-------
!NEUTRAL WIND SOURCE TERMS FOR POTENTIAL EQUATION, SIMILAR TO ABOVE BLOCK OF CODE
J1=0d0    !so this div is only perp components
if (flagswap==1) then
  J2=-1*sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)     !wind x2 current, note that all workers already have a copy of this.  
  J3=sigH*vn3*B1(1:lx1,1:lx2,1:lx3)+sigP*vn2*B1(1:lx1,1:lx2,1:lx3)     !wind x3 current
else
  J2=sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)     !wind x2 current
  J3=sigH*vn3*B1(1:lx1,1:lx2,1:lx3)-sigP*vn2*B1(1:lx1,1:lx2,1:lx3)     !wind x3 current
end if
!    divJperp=div3D(J1,J2,J3,x,1,lx1,1,lx2,1,lx3)
J1halo(1:lx1,1:lx2,1:lx3)=J1
J2halo(1:lx1,1:lx2,1:lx3)=J2
J3halo(1:lx1,1:lx2,1:lx3)=J3

J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)

J2halo(0,1:lx2,1:lx3)=J2halo(1,1:lx2,1:lx3)
J2halo(lx1+1,1:lx2,1:lx3)=J2halo(lx1,1:lx2,1:lx3)

J3halo(0,1:lx2,1:lx3)=J3halo(1,1:lx2,1:lx3)
J3halo(lx1+1,1:lx2,1:lx3)=J3halo(lx1,1:lx2,1:lx3)

call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 
call halo(J2halo,1,tagJ2,x%flagper)
call halo(J3halo,1,tagJ3,x%flagper)

!ZZZ - CHECK FOR AND SET GLOBAL BOUNDARIES
if (iddown==-1) then
  J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  J2halo(1:lx1,0,1:lx3)=J2halo(1:lx1,1,1:lx3)
  J3halo(1:lx1,0,1:lx3)=J3halo(1:lx1,1,1:lx3)
end if
if (idup==lid2) then
  J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
  J2halo(1:lx1,lx2+1,1:lx3)=J2halo(1:lx1,lx2,1:lx3)
  J3halo(1:lx1,lx2+1,1:lx3)=J3halo(1:lx1,lx2,1:lx3)
end if
if (.not. x%flagper) then
  if (idleft==-1) then
    J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
    J2halo(1:lx1,1:lx2,0)=J2halo(1:lx1,1:lx2,1)
    J3halo(1:lx1,1:lx2,0)=J3halo(1:lx1,1:lx2,1)
  end if
  if (idright==lid3) then
    J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
    J2halo(1:lx1,1:lx2,lx3+1)=J2halo(1:lx1,1:lx2,lx3)
    J3halo(1:lx1,1:lx2,lx3+1)=J3halo(1:lx1,1:lx2,lx3)
  end if
end if

divtmp=div3D(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),J2halo(0:lx1+1,0:lx2+1,0:lx3+1), &
             J3halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
srcterm=srcterm+divtmp(1:lx1,1:lx2,1:lx3)
!-------


!    !ZZZ - DEBUG BY GETTING THE ENTIRE SOURCETERM ARRAY
!    call gather_send(srcterm,tagsrc)


!!!!!!!!
!-----AT THIS POINT WE MUST DECIDE WHETHER TO DO AN INTEGRATED SOLVE OR A 2D FIELD-RESOLVED SOLVED
!-----DECIDE BASED ON THE SIZE OF THE X2 DIMENSION
if (lx2/=1) then    !either field-resolved 3D or integrated 2D solve for 3D domain
  if (potsolve == 1) then    !2D, field-integrated solve
    !-------
    !INTEGRATE CONDUCTANCES AND CAPACITANCES FOR SOLVER COEFFICIENTS 
    integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)
    sigintegral=integral3D1(integrand,x,1,lx1)    !no haloing required for a field-line integration
    SigPint2=sigintegral(lx1,:,:)

    integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)
    sigintegral=integral3D1(integrand,x,1,lx1)
    SigPint3=sigintegral(lx1,:,:)

    integrand=x%h1(1:lx1,1:lx2,1:lx3)*sigH
    sigintegral=integral3D1(integrand,x,1,lx1)
    SigHint=sigintegral(lx1,:,:)
    
    sigintegral=integral3D1(incap,x,1,lx1)
    incapint=sigintegral(lx1,:,:)
    !-------


    !PRODUCE A FIELD-INTEGRATED SOURCE TERM
    if (flagdirich /= 1) then   !Neumann conditions; incorporate a source term and execute the solve
      !-------
      integrand=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)*srcterm
      sigintegral=integral3D1(integrand,x,1,lx1)
      srctermint=sigintegral(lx1,:,:)
      srctermint=srctermint+x%h2(lx1,1:lx2,1:lx3)*x%h3(lx1,1:lx2,1:lx3)*Vmaxx1slab- &
                                 x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab     !workers don't have access to boundary conditions, unless root sends
      !-------


      !RADD--- ROOT NEEDS TO PICK UP *INTEGRATED* SOURCE TERMS AND COEFFICIENTS FROM WORKERS
      call gather_send(srctermint,tagsrc)
      call gather_send(incapint,tagincapint)
      call gather_send(SigPint2,tagSigPint2)
      call gather_send(SigPint3,tagSigPint3)
      call gather_send(SigHint,tagSigHint)
      v2=vs2(1:lx1,1:lx2,1:lx3,1); v3=vs3(1:lx1,1:lx2,1:lx3,1);
      v2slab=v2(lx1,:,:); v3slab=v3(lx1,:,:)
      call gather_send(v2slab,tagv2electro)
      call gather_send(v3slab,tagv3electro)
!          v2slab=vs2(lx1,1:lx2,1:lx3,1); v3slab=vs3(lx1,1:lx2,1:lx3,1);    !need to pick out the ExB drift here (i.e. the drifts from highest altitudes); but this is only valid for Cartesian, so it's okay for the foreseeable future


      call elliptic_workers()    !workers do not need any specific info about      proglem.  

    else     !Dirichlet conditions - since this is field integrated we just copy BCs specified by user to other locations along field line

    end if

!
  else    !resolved 3D solve
    !ZZZ - conductivities need to be properly scaled here...  So does the source term...  Maybe leave as broken for now since I don't really plan to use this code


    !RADD--- ROOT NEEDS TO PICK UP FIELD-RESOLVED SOURCE TERM AND COEFFICIENTS FROM WORKERS
    call gather_send(sigPscaled,tagsigP)
    call gather_send(sigHscaled,tagsigH)
    call gather_send(sig0scaled,tagsig0)
    call gather_send(srcterm,tagsrc)

    call elliptic_workers()

  end if
else   !lx1=1 so do a field-resolved 2D solve over x1,x3


  !-------
  !PRODUCE SCALED CONDUCTIVITIES TO PASS TO SOLVER, ALSO SCALED SOURCE TERM
  sig0scaled=x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)*sig0
  if (flagswap==1) then
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)*sigP    !remember to swap 2-->3
  else
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
  end if
  srcterm=srcterm*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
  !-------

  !RADD--- NEED TO GET THE RESOLVED SOURCE TERMS AND COEFFICIENTS FROM WORKERS
  call gather_send(sigPscaled,tagsigP)
  call gather_send(sig0scaled,tagsig0)
  call gather_send(srcterm,tagsrc)

  call elliptic_workers()
end if
!    print *, 'MUMPS time:  ',tfin-tstart
!!!!!!!!!


!RADD--- ROOT NEEDS TO PUSH THE POTENTIAL BACK TO ALL WORKERS FOR FURTHER PROCESSING (BELOW)
call bcast_recv(Phi,tagPhi)


!-------
!STORE PREVIOUS TIME TOTAL FIELDS BEFORE UPDATING THE ELECTRIC FIELDS WITH NEW POTENTIAL (OLD FIELDS USED TO CALCULATE POLARIZATION CURRENT)
E1prev=E1
E2prev=E2
E3prev=E3
!-------


!-------
!CALCULATE PERP FIELDS FROM POTENTIAL
!      E20all=grad3D2(-1d0*Phi0all,dx2(1:lx2))    !causes major memory leak. maybe from arithmetic statement argument? Left here as a 'lesson learned' (or is it a gfortran bug...)
!      E30all=grad3D3(-1d0*Phi0all,dx3all(1:lx3all))
Phi=-1d0*Phi
E2=grad3D2(Phi,x,1,lx1,1,lx2,1,lx3)    !no haloing required
!    E3=grad3D3(Phi,x,1,lx1,1,lx2,1,lx3)    !needs to be haloed

J1halo(1:lx1,1:lx2,1:lx3)=Phi
J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)
call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 

if (iddown==-1) then
  J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
end if
if (idup==lid2) then
  J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
end if
if (.not. x%flagper) then
  if (idleft==-1) then
    J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
  end if
  if (idright==lid3) then
    J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
  end if
end if
  
divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
E3=divtmp(1:lx1,1:lx2,1:lx3)
Phi=-1d0*Phi   !put things back for later use
!--------


!    !R-------
!    !JUST TO JUDGE THE IMPACT OF MI COUPLING
!    print *, 'Max integrated inertial capacitance:  ',maxval(incapint)
!    print *, 'Max integrated Pedersen conductance (includes metric factors):  ',maxval(SigPint2)
!    print *, 'Max integrated Hall conductance (includes metric factors):  ',minval(SigHint), maxval(SigHint)
!!    print *, 'Max E2,3 BG and response values are:  ',maxval(abs(E02)), maxval(abs(E03)), maxval(abs(E2)),maxval(abs(E3))
!    print *, 'Max E2,3 BG and response values are:  ',maxval(E02), maxval(E03),maxval(E2),maxval(E3)
!    print *, 'Min E2,3 BG and response values are:  ',minval(E02), minval(E03),minval(E2),minval(E3)
!    !R-------


!--------
!ADD IN BACKGROUND FIELDS BEFORE DISTRIBUTING TO WORKERS
E2=E2+E02
E3=E3+E03
!--------


!--------
!COMPUTE TIME DERIVATIVE NEEDED FOR POLARIZATION CURRENT.  ONLY DO THIS IF WE HAVE SPECIFIC NONZERO INERTIAL CAPACITANCE
if (maxval(incap) > 1d-6) then
  grad2E=grad3D2(E2,x,1,lx1,1,lx2,1,lx3)
!      grad3E=grad3D3(E2,x,1,lx1,1,lx2,1,lx3)

  J1halo(1:lx1,1:lx2,1:lx3)=E2
  J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
  J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)
  call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 

  if (iddown==-1) then
    J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  end if
  if (idup==lid2) then
    J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
  end if
  if (.not. x%flagper) then
    if (idleft==-1) then
      J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
    end if
    if (idright==lid3) then
      J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
    end if
  end if

  divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad3E=divtmp(1:lx1,1:lx2,1:lx3)

  DE2Dt=(E2-E2prev)/dt+v2*grad2E+v3*grad3E
  grad2E=grad3D2(E3,x,1,lx1,1,lx2,1,lx3)
!      grad3E=grad3D3(E3,x,1,lx1,1,lx2,1,lx3)

  J1halo(1:lx1,1:lx2,1:lx3)=E3
  J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
  J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)
  call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 

  if (iddown==-1) then
    J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
  end if
  if (idup==lid2) then
    J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
  end if
  if (.not. x%flagper) then
    if (idleft==-1) then
      J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
    end if
    if (idright==lid3) then
      J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
    end if
  end if

  divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad3E=divtmp(1:lx1,1:lx2,1:lx3)

  DE3Dt=(E3-E3prev)/dt+v2*grad2E+v3*grad3E
  J1pol=0d0
  J2pol=incap*DE2Dt
  J3pol=incap*DE3Dt
else       !pure electrostatic solve was done
  DE2Dt=0d0
  DE3Dt=0d0
  J1pol=0d0 
  J2pol=0d0
  J3pol=0d0      
end if
!--------

!-------
if (flagswap==1) then
  J2=sigP*E2+sigH*E3    !BG field already added to E above
  J3=-1*sigH*E2+sigP*E3
else
  J2=sigP*E2-sigH*E3    !BG field already added to E above
  J3=sigH*E2+sigP*E3
end if

!WHAT I THINK THE NEUTRAL WIND CURRENTS SHOULD BE IN 2D
if (flagswap==1) then
  J2=J2-sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
  J3=J3+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)+sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
else
  J2=J2+sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
  J3=J3+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)-sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
end if
!-------


!!!!!!!!
!NOW DEAL WITH THE PARALLEL FIELDS AND ALL CURRENTS
if (lx2/=1 .and. potsolve ==1) then    !we did a field-integrated solve above

  !-------
  !NOTE THAT A DIRECT E1ALL CALCULATION WILL GIVE ZERO, SO USE INDIRECT METHOD, AS FOLLOWS
  J1=0d0    !a placeholder so that only the perp divergence is calculated - will get overwritten later.
!      divJperp=div3D(J1,J2,J3,x,1,lx1,1,lx2,1,lx3)
  J1halo(1:lx1,1:lx2,1:lx3)=J1
  J2halo(1:lx1,1:lx2,1:lx3)=J2
  J3halo(1:lx1,1:lx2,1:lx3)=J3

  J1halo(0,1:lx2,1:lx3)=J1halo(1,1:lx2,1:lx3)
  J1halo(lx1+1,1:lx2,1:lx3)=J1halo(lx1,1:lx2,1:lx3)

  J2halo(0,1:lx2,1:lx3)=J2halo(1,1:lx2,1:lx3)
  J2halo(lx1+1,1:lx2,1:lx3)=J2halo(lx1,1:lx2,1:lx3)

  J3halo(0,1:lx2,1:lx3)=J3halo(1,1:lx2,1:lx3)
  J3halo(lx1+1,1:lx2,1:lx3)=J3halo(lx1,1:lx2,1:lx3)

  call halo(J1halo,1,tagJ1,x%flagper)    !I'm kind of afraid to only halo a single point... 
  call halo(J2halo,1,tagJ2,x%flagper)
  call halo(J3halo,1,tagJ3,x%flagper)

  !ZZZ - SET GLOBAL BOUNDARY HERE
  if (iddown==-1) then
    J1halo(1:lx1,0,1:lx3)=J1halo(1:lx1,1,1:lx3)
    J2halo(1:lx1,0,1:lx3)=J2halo(1:lx1,1,1:lx3)
    J3halo(1:lx1,0,1:lx3)=J3halo(1:lx1,1,1:lx3)
  end if
  if (idup==lid2) then
    J1halo(1:lx1,lx2+1,1:lx3)=J1halo(1:lx1,lx2,1:lx3)
    J2halo(1:lx1,lx2+1,1:lx3)=J2halo(1:lx1,lx2,1:lx3)
    J3halo(1:lx1,lx2+1,1:lx3)=J3halo(1:lx1,lx2,1:lx3)
  end if
  if (.not. x%flagper) then
    if (idleft==-1) then
      J1halo(1:lx1,1:lx2,0)=J1halo(1:lx1,1:lx2,1)
      J2halo(1:lx1,1:lx2,0)=J2halo(1:lx1,1:lx2,1)
      J3halo(1:lx1,1:lx2,0)=J3halo(1:lx1,1:lx2,1)
    end if
    if (idright==lid3) then
      J1halo(1:lx1,1:lx2,lx3+1)=J1halo(1:lx1,1:lx2,lx3)
      J2halo(1:lx1,1:lx2,lx3+1)=J2halo(1:lx1,1:lx2,lx3)
      J3halo(1:lx1,1:lx2,lx3+1)=J3halo(1:lx1,1:lx2,lx3)
    end if
  end if

  divtmp=div3D(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),J2halo(0:lx1+1,0:lx2+1,0:lx3+1), &
               J3halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  divJperp=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)*divtmp(1:lx1,1:lx2,1:lx3)
  if (flagdirich /= 1) then     !Neumann conditions, this is boundary location-agnostic since both bottom and top FACs are known - they have to  be loaded into VVmaxx1 and Vminx1.  For numerical purposes we prefer to integrate from the location of nonzero current (usually highest altitude in open grid).  
    if (gridflag==0) then     !closed dipole grid, really would be best off integrating from the source hemisphere
!          print *,  'Closed dipole grid; integration starting at max x1...', minval(Vmaxx1slab), &
!                         maxval(Vmaxx1slab)
      if (sourcemlat>=0d0) then    !integrate from northern hemisphere
!            print *, 'Source is in northern hemisphere (or there is no source)...'
        J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,:,:))
        end do
      else
!            print *, 'Source in southern hemisphere...'
        J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,:,:))
        end do
      end if
    elseif (gridflag==1) then    !this would be an inverted grid, this max altitude corresponds to the min value of x1
!          print *,  'Inverted grid; integration starting at min x1...',minval(Vminx1slab), maxval(Vminx1slab)
      J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
      do ix1=1,lx1  
        J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                         (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,:,:))
      end do
    else        !minx1 is at teh bottom of the grid to integrate from max x1
!          print *,  'Non-inverted grid; integration starting at max x1...', minval(Vmaxx1slab),  maxval(Vmaxx1slab)
      J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
      do ix1=1,lx1
        J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                         (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,:,:))
      end do
    end if
!        if (gridflag==2) then    !for a cartesian grid in the northern hemisphere (assumed) we have the x1-direction being against the magnetic field...
!          J1=-1d0*J1     !ZZZ - very questionable
!        end if
  else      !Dirichlet conditions - we need to integrate from the ***lowest altitude*** (where FAC is known to be zero, note this is not necessarilty the logical bottom of the grid), upwards (to where it isn't)
    if (gridflag/=2) then    !inverted grid (logical top is the lowest altitude)
!          print *, 'Inverted grid detected - integrating logical top downward to compute FAC...'
      J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current
      do ix1=1,lx1
        J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                         (J1(ix1,:,:))    !FAC AT TOP ASSUMED TO BE ZERO
      end do
    else      !non-inverted grid (logical bottom is the lowest altitude - so integrate normy)
!          print *, 'Non-inverted grid detected - integrating logical bottom to top to compute FAC...'
      J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
      do ix1=1,lx1
        J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                         (-1d0*J1(ix1,:,:))    !FAC AT THE BOTTOM ASSUMED TO BE ZERO
      end do
    end if
  end if
  E1=J1/sig0
  !-------

else   !we resolved the field line (either 2D solve or full 3D) so just differentiate normally

  !-------
  Phi=-1d0*Phi
  E1=grad3D1(Phi,x,1,lx1,1,lx2,1,lx3)    !no haloing required since x1-derivative
  Phi=-1d0*Phi
  J1=sig0*E1
  !-------

end if
!!!!!!!!!


!    !R-------
!    print *, 'Max topside FAC (abs. val.) computed to be:  ',maxval(abs(J1(1,:,:)))    !ZZZ - this rey needsz to be current at the "top"
!    print *, 'Max polarization J2,3 (abs. val.) computed to be:  ',maxval(abs(J2pol)),  maxval(abs(J3pol))
!    print *, 'Max conduction J2,3  computed to be:  ',maxval(J2), maxval(J3)
!    print *, 'Min conduction J2,3  computed to be:  ',minval(J2),  minval(J3)
!    print *, 'Max conduction J1 (abs. val.) computed to be:  ',maxval(abs(J1))
!    print *, 'flagswap:  ',flagswap
!    !R-------


!-------
!GRAND TOTAL FOR THE CURRENT DENSITY:  TOSS IN POLARIZATION CURRENT SO THAT OUTPUT FILES ARE CONSISTENT
J1=J1+J1pol
J2=J2+J2pol
J3=J3+J3pol
!-------

end subroutine potential_workers_mpi

end module potential_comm

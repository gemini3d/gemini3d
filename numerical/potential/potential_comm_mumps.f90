module potential_comm

!! THIS MODULE IS MEANT TO WORK WITH THE MUMPS 2D
!! INTEGRATED SOLVER IF THE GRID IS 3D, OR A FIELD-RESOLVED
!! SOLVER IF THE GRID IS 2D (MUMPS CAN'T HANDLE 3D VERY WELL).
!!
!! NOTE THAT ONLY THE CURVILINEAR FUNCTION ARE UP TO DATE.

use mpi, only: mpi_integer, mpi_comm_world, mpi_status_ignore

use phys_consts, only: wp, pi, lsp, debug
use grid, only: flagswap, gridflag
use mesh, only: curvmesh
use collisions, only: conductivities, capacitance
use calculus, only: div3d, integral3d1, grad3d1, grad3d2, grad3d3, integral3d1_curv_alt
use potentialBCs_mumps, only: potentialbcs2D, potentialbcs2D_fileinput
use potential_mumps, only: potential3D_fieldresolved_decimate, &
                            potential2D_fieldresolved, &
                            potential2D_polarization, &
                            potential2D_polarization_periodic
use PDEelliptic, only: elliptic_workers

use mpimod, only: lid, lid2, lid3, myid, myid2, myid3, tage01, tage02, tage03, tagflagdirich, tagincapint, &
  tagsrc, tagv2electro, tagv3electro, tagvmaxx1, tagvminx1, tagj1, tagj2, tagj3, tagphi, tagsig0, &
  tagsigh, tagsighint, tagsigp, tagsigpint2, tagsigpint3, &
  bcast_send, bcast_recv, gather_recv, gather_send, halo

implicit none
private
public :: electrodynamics, halo_pot


!! OVERLOAD THE SOLVERS FOR DEALING WITH THE CURVILINEAR MESHES
!! NOTE WORKER SUBROUTINE DOES NOT NEED TO BE CHANGED/OVERLOADED
interface electrodynamics
  module procedure electrodynamics_curv
end interface electrodynamics

interface potential_root_mpi
  module procedure potential_root_mpi_curv
end interface potential_root_mpi

interface ! potential_worker.f90
module subroutine potential_workers_mpi(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,vn2,vn3,sourcemlat,B1,x, &
  potsolve,flagcap, &
  E1,E2,E3,J1,J2,J3)

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
integer, intent(in) :: flagcap

real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
end subroutine potential_workers_mpi
end interface

interface ! potential_root.f90
module subroutine potential_root_mpi_curv(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,vn2,vn3,sourcemlat,B1,x, &
  potsolve,flagcap, &
  E1,E2,E3,J1,J2,J3, &
  Phiall,flagE0file,dtE0,E0dir,ymd,UTsec)

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
integer, intent(in) :: flagcap

real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
real(wp), dimension(:,:,:), intent(inout) :: Phiall   !not good form, but I'm lazy...  Forgot what I meant by this...

integer, intent(in) :: flagE0file
real(wp), intent(in) :: dtE0
character(*), intent(in) :: E0dir
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
end subroutine potential_root_mpi_curv
end interface

contains


subroutine electrodynamics_curv(it,t,dt,nn,vn2,vn3,Tn,sourcemlat,ns,Ts,vs1,B1,vs2,vs3,x, &
                         potsolve,flagcap, &
                         E1,E2,E3,J1,J2,J3,Phiall, &
                         flagE0file,dtE0,E0dir,ymd,UTsec)

!! THIS IS A WRAPPER FUNCTION FOR THE ELECTRODYANMICS
!! PART OF THE MODEL.  BOTH THE ROOT AND WORKER PROCESSES
!! CALL THIS SAME SUBROUTINE, WHEN THEN BRANCHES INTO
!! DIFFERENT TASKS FOR EACH AFTER ALL COMPUTE CONDUCTIVITIES
!! AND INERTIAL CAPACITANCE.
!!
!! NOTE THAT THE ALLOCATION STATUS
!! OF THE ALL VARIABLES FOR THE WORKERS WILL BE UNALLOCATED.
!! THIS code requires FORTRAN >= 2003 STANDARD.

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
real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall
!! inout since it may not be allocated or deallocated in this procedure

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
    if (debug) print *, 'Conductivities and capacitance for time step:  ',t,' took ',tfin-tstart,' seconds...'
  else
    if (debug) print *, 'Conductivities for time step:  ',t,' took ',tfin-tstart,' seconds...'
  end if
end if

if (potsolve == 1 .or. potsolve == 3) then    !electrostatic solve or electrostatic alt. solve
  call cpu_time(tstart)

  if (myid/=0) then    !role-specific communication pattern (all-to-root-to-all), workers initiate with sends
     call potential_workers_mpi(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,vn2,vn3,sourcemlat,B1,x, &
                            potsolve,flagcap, &
                            E1,E2,E3,J1,J2,J3)
  else
    call potential_root_mpi(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,vn2,vn3,sourcemlat,B1,x, &
                              potsolve,flagcap, &
                              E1,E2,E3,J1,J2,J3, &
                              Phiall,flagE0file,dtE0,E0dir,ymd,UTsec)
  end if

  !DRIFTS - NEED TO INCLUDE ELECTRIC, WIND-DRIVEN, AND GRAVITATIONAL???
!  if (lx2/=1) then    !full 3D solve, go with the regular formulas
  if(flagswap/=1) then
    do isp=1,lsp
      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+muPvn(:,:,:,isp)*vn2-muHvn(:,:,:,isp)*vn3
      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+muHvn(:,:,:,isp)*vn2+muPvn(:,:,:,isp)*vn3
    end do
  else                !flip signs on the cross products in 2D.  Note that due to dimension shuffling E2,3 mapping is already handled
    do isp=1,lsp
      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2+muH(:,:,:,isp)*E3+muPvn(:,:,:,isp)*vn2+muHvn(:,:,:,isp)*vn3
      vs3(1:lx1,1:lx2,1:lx3,isp)=-muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3-muHvn(:,:,:,isp)*vn2+muPvn(:,:,:,isp)*vn3
    end do
  end if

!    do isp=1,lsp
       !! To leading order the ion drifts do not include the polarization parts,
       !! otherwise it may mess up polarization convective term in the electrodynamics solver...
!      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE2Dt
!      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE3Dt
!    end do

  call cpu_time(tfin)

  if (myid==0) then
    if (debug) print *, 'Potential solution for time step:  ',t,' took ',tfin-tstart,' seconds...'
    if (debug) print *, 'Min and max root drift values:  ',minval(vs2),maxval(vs2), minval(vs3),maxval(vs3)
  end if

else if (potsolve == 2) then  !inductive form of model, could this be subcycled to speed things up?
  !Do nothing for now...
else   !null solve; just force everything to zero
  E1=0d0; E2=0d0; E3=0d0; J1=0d0; J2=0d0; J3=0d0;
  vs2=0d0; vs3=0d0;
end if

end subroutine electrodynamics_curv


subroutine halo_pot(parmhalo,tagcurrent,flagper,flagdegrade)

!THIS SUBROUTINE REPLICATES A COMMON MESSAGE PASSING SCHEME USED IN THE COMPUTATION
!OF ELECTRODYNAMICS PARAMETERS THAT RESULT FROM DERIVATIVE (WHICH REQUIRE HALOING)

real(wp), intent(inout), dimension(-1:,-1:,-1:) :: parmhalo
integer, intent(in) :: tagcurrent
logical, intent(in) :: flagper
logical, intent(in) :: flagdegrade    !whether or not to degrade edge derivatives to first order in x2

integer :: lx1,lx2,lx3
integer :: idleft,idright,iddown,idup


lx1=size(parmhalo,1)-4
lx2=size(parmhalo,2)-4
lx3=size(parmhalo,3)-4

idleft=myid3-1
idright=myid3+1
iddown=myid2-1
idup=myid2+1

parmhalo(0,1:lx2,1:lx3)=parmhalo(1,1:lx2,1:lx3)
parmhalo(lx1+1,1:lx2,1:lx3)=parmhalo(lx1,1:lx2,1:lx3)

call halo(parmhalo,1,tagcurrent,flagper)     !this particular type of message passing only needs a single ghost cell

if (iddown==-1) then
  if (flagdegrade .and. lx2>1) then     !for whatever reason this fails ctest without checking lx2>1
    parmhalo(1:lx1,0,1:lx3)=-1*parmhalo(1:lx1,2,1:lx3)+2*parmhalo(1:lx1,1,1:lx3)
  else
    parmhalo(1:lx1,0,1:lx3)=parmhalo(1:lx1,1,1:lx3)
  end if
end if
if (idup==lid2) then
  if (flagdegrade .and. lx2>1) then
    parmhalo(1:lx1,lx2+1,1:lx3)=2*parmhalo(1:lx1,lx2,1:lx3)-parmhalo(1:lx1,lx2-1,1:lx3)
  else
    parmhalo(1:lx1,lx2+1,1:lx3)=parmhalo(1:lx1,lx2,1:lx3)
  end if
end if
if (.not. flagper) then              !musn't overwrite ghost cells if perioidc is chosen
  if (idleft==-1) then
    parmhalo(1:lx1,1:lx2,0)=parmhalo(1:lx1,1:lx2,1)
  end if
  if (idright==lid3) then
    parmhalo(1:lx1,1:lx2,lx3+1)=parmhalo(1:lx1,1:lx2,lx3)
  end if
end if

end subroutine halo_pot


end module potential_comm

module potential_comm

!! THIS MODULE IS MEANT TO WORK WITH THE MUMPS 2D
!! INTEGRATED SOLVER IF THE GRID IS 3D, OR A FIELD-RESOLVED
!! SOLVER IF THE GRID IS 2D (MUMPS CAN'T HANDLE 3D VERY WELL).
!!
!! NOTE THAT ONLY THE CURVILINEAR FUNCTION ARE UP TO DATE.

use, intrinsic :: ieee_arithmetic

use phys_consts, only: wp, pi, lsp, debug, ms, qs, kB
use grid, only: flagswap, gridflag, lx1,lx2all,lx3all, g1,g2,g3
use mesh, only: curvmesh
use collisions, only: conductivities, capacitance
use calculus, only: div3d, integral3d1, grad3d1, grad3d2, grad3d3, integral3d1_curv_alt
use potentialBCs_mumps, only: potentialbcs2D, potentialbcs2D_fileinput, compute_rootBGEfields
use potential_mumps, only: potential3D_fieldresolved_decimate, &
                            potential2D_fieldresolved, &
                            potential2D_polarization, &
                            potential2D_polarization_periodic
use PDEelliptic, only: elliptic_workers
use mpimod, only: mpi_integer, mpi_comm_world, mpi_status_ignore, &
mpi_cfg, tag=>gemini_mpi, &
bcast_send, bcast_recv, gather_recv, gather_send, halo
use config, only: gemini_cfg

implicit none (type, external)
private
public :: electrodynamics, halo_pot, potential_sourceterms, pot2perpfield, velocities, get_BGEfields, &
            acc_perpconductioncurrents,acc_perpwindcurrents,acc_perpgravcurrents

external :: mpi_send, mpi_recv

!! overloading to deal with vestigial cartesian->curvilinear code
interface electrodynamics
  module procedure electrodynamics_curv
end interface electrodynamics

interface potential_root_mpi
  module procedure potential_root_mpi_curv
end interface potential_root_mpi

interface ! potential_worker.f90
  module subroutine potential_workers_mpi(it,t,dt,sig0,sigP,sigH,sigPgrav,sigHgrav, &
                                            muP,muH, &
                                            incap,vs2,vs3,vn2,vn3,cfg,B1,ns,Ts,x, &
                                            E1,E2,E3,J1,J2,J3)

  integer, intent(in) :: it
  real(wp), intent(in) :: t,dt
  real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH,sigPgrav,sigHgrav
  real(wp), dimension(:,:,:,:), intent(in) :: muP,muH
  real(wp), dimension(:,:,:), intent(in) ::  incap
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3
  real(wp), dimension(:,:,:), intent(in) ::  vn2,vn3
  type(gemini_cfg), intent(in) :: cfg
  real(wp), dimension(-1:,-1:,-1:), intent(in) ::  B1
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  ns,Ts

  type(curvmesh), intent(in) :: x

  real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
  end subroutine potential_workers_mpi
end interface

interface ! potential_root.f90
  module subroutine potential_root_mpi_curv(it,t,dt,sig0,sigP,sigH,sigPgrav,sigHgrav, &
                                              muP,muH, &
                                              incap,vs2,vs3,vn2,vn3,cfg,B1,ns,Ts,x, &
                                              E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)

  integer, intent(in) :: it
  real(wp), intent(in) :: t,dt
  real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH,sigPgrav,sigHgrav
  real(wp), dimension(:,:,:,:), intent(in) :: muP,muH
  real(wp), dimension(:,:,:), intent(in) ::  incap
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3
  real(wp), dimension(:,:,:), intent(in) ::  vn2,vn3
  type(gemini_cfg), intent(in) :: cfg
  real(wp), dimension(-1:,-1:,-1:), intent(in) ::  B1
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  ns,Ts

  type(curvmesh), intent(in) :: x

  real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
  real(wp), dimension(:,:,:), intent(inout) :: Phiall   !not good form, but I'm lazy...  Forgot what I meant by this...
  integer, dimension(3), intent(in) :: ymd
  real(wp), intent(in) :: UTsec
  end subroutine potential_root_mpi_curv
end interface

contains


subroutine electrodynamics_curv(it,t,dt,nn,vn2,vn3,Tn,cfg,ns,Ts,vs1,B1,vs2,vs3,x, &
                         E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)

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
type(gemini_cfg), intent(in) :: cfg
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  vs2,vs3
type(curvmesh), intent(in) :: x
real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall
!! inout since it may not be allocated or deallocated in this procedure
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec

real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: sig0,sigP,sigH,sigPgrav,sigHgrav
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,1:size(ns,4)) :: muP,muH,nusn
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
call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)


! Error checking for cap. vs. grid types
if (cfg%flagcap/=0) then
  call capacitance(ns,B1,cfg,incap)    !> full cfg needed for optional inputs...
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
if (mpi_cfg%myid==0) then
  if (cfg%flagcap/=0) then
    if (debug) print *, 'Conductivities and capacitance for time step:  ',t,' took ',tfin-tstart,' seconds...'
  else
    if (debug) print *, 'Conductivities for time step:  ',t,' took ',tfin-tstart,' seconds...'
  end if
end if


if (cfg%potsolve == 1 .or. cfg%potsolve == 3) then    !electrostatic solve or electrostatic alt. solve
  call cpu_time(tstart)

  ! Execute solution for ionospheric potential
  if (mpi_cfg%myid/=0) then    
    !! role-specific communication pattern (all-to-root-to-all), workers initiate with sends
     call potential_workers_mpi(it,t,dt,sig0,sigP,sigH,sigPgrav,sigHgrav,muP,muH,incap,vs2,vs3, &
                                 vn2,vn3,cfg,B1,ns,Ts,x,E1,E2,E3,J1,J2,J3)
  else
    call potential_root_mpi(it,t,dt,sig0,sigP,sigH,sigPgrav,sigHgrav,muP,muH,incap,vs2,vs3,vn2,vn3,cfg,B1,ns,Ts,x, &
                              E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)
  end if

  ! Compute velocity from mobilities and fields
  call velocities(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)

  call cpu_time(tfin)

  if (mpi_cfg%myid==0) then
    if (debug) print *, 'Potential solution for time step:  ',t,' took ',tfin-tstart,' seconds...'
    if (debug) print *, 'Min and max root drift values:  ',minval(vs2),maxval(vs2), minval(vs3),maxval(vs3)
  end if

else if (cfg%potsolve == 2) then  !inductive form of model, could this be subcycled to speed things up?
  error stop 'Inductive solves are not supported yet...'
else   !null solve; just force everything to zero
  E1=0d0; E2=0d0; E3=0d0; J1=0d0; J2=0d0; J3=0d0;
  vs2=0d0; vs3=0d0;
end if

end subroutine electrodynamics_curv


subroutine velocities(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,flaggravdrift,flagdiamagnetic,vs2,vs3)

!> compute steady state drifts resulting from a range of forces.  Can be used
!   by both root and worker processes

real(wp), dimension(:,:,:,:), intent(in) :: muP,muH,nusn
real(wp), dimension(:,:,:), intent(in) :: E2,E3,vn2,vn3
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts      !> these must have ghost cells
type(curvmesh), intent(in) :: x
logical, intent(in) :: flaggravdrift
logical, intent(in) :: flagdiamagnetic
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: vs2,vs3   !> these have ghost cells

integer :: lx1,lx2,lx3,lsp,isp
real(wp), dimension(-1:size(E2,1)+2,-1:size(E2,2)+2,-1:size(E2,3)+2) :: pressure
real(wp), dimension(0:size(E2,1)+1,0:size(E2,2)+1,0:size(E2,3)+1) :: gradpx,gradpy

!> sizes from the mobility coefficients
lx1=size(muP,1)
lx2=size(muP,2)
lx3=size(muP,3)
lsp=size(muP,4)


!! FIXME:  Do we really need separate wind mobility or can we just compute off electrical mobility as done with gravity.  This is necessary because we are not storing the collision frequencies...  For adding pressure terms we may want to just use collision freq. to avoid a bunch of slightly different mobility arrays.
!> electric field and wind terms for ion drifts
if(flagswap/=1) then
  do isp=1,lsp
    vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+ &
                      (muP(:,:,:,isp)*vn2-muH(:,:,:,isp)*vn3)*(ms(isp)*nusn(:,:,:,isp)/qs(isp))
    vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ &
                      (muH(:,:,:,isp)*vn2+muP(:,:,:,isp)*vn3)*ms(isp)*nusn(:,:,:,isp)/qs(isp)
  end do
else                !flip signs on the cross products in 2D.  Note that due to dimension shuffling E2,3 mapping is already handled.
  do isp=1,lsp
    vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2+muH(:,:,:,isp)*E3+ &
                      (muP(:,:,:,isp)*vn2+muH(:,:,:,isp)*vn3)*ms(isp)*nusn(:,:,:,isp)/qs(isp)
    vs3(1:lx1,1:lx2,1:lx3,isp)=-muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ &
                      (-1._wp*muH(:,:,:,isp)*vn2+muP(:,:,:,isp)*vn3)*ms(isp)*nusn(:,:,:,isp)/qs(isp)
  end do
end if

!> Pressure/diamagnetic terms (if required)
if (flagdiamagnetic) then
  if(flagswap/=1) then
    do isp=1,lsp
      pressure=ns(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)    ! compute pressure from n,T
      call halo_pot(pressure,tag%pressure,x%flagper,.false.)             ! boundary fill via haloing
      gradpx=grad3D2(pressure,x,0,lx1+1,0,lx2+1,0,lx3+1)                 ! compute gradient x2,x3 components
      gradpy=grad3D3(pressure,x,0,lx1+1,0,lx2+1,0,lx3+1)
      vs2(1:lx1,1:lx2,1:lx3,isp)=vs2(1:lx1,1:lx2,1:lx3,isp) &
                -muP(1:lx1,1:lx2,1:lx3,isp)/ns(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradpx(1:lx1,1:lx2,1:lx3) &
                +muH(1:lx1,1:lx2,1:lx3,isp)/ns(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradpy(1:lx1,1:lx2,1:lx3)
      vs3(1:lx1,1:lx2,1:lx3,isp)=vs3(1:lx1,1:lx2,1:lx3,isp) &
                -muH(1:lx1,1:lx2,1:lx3,isp)/ns(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradpx(1:lx1,1:lx2,1:lx3) &
                -muP(1:lx1,1:lx2,1:lx3,isp)/ns(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradpy(1:lx1,1:lx2,1:lx3)
    end do
  else             !coordinates have been swapped so flip signs on the cross product (Hall) terms
    do isp=1,lsp
      pressure=ns(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)    ! compute pressure from n,T
      call halo_pot(pressure,tag%pressure,x%flagper,.false.)             ! boundary fill via haloing
      gradpx=grad3D2(pressure,x,0,lx1+1,0,lx2+1,0,lx3+1)                 ! compute gradient x2,x3 components
      gradpy=grad3D3(pressure,x,0,lx1+1,0,lx2+1,0,lx3+1)
      vs2(1:lx1,1:lx2,1:lx3,isp)=vs2(1:lx1,1:lx2,1:lx3,isp) &
                -muP(1:lx1,1:lx2,1:lx3,isp)/ns(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradpx(1:lx1,1:lx2,1:lx3) &
                -muH(1:lx1,1:lx2,1:lx3,isp)/ns(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradpy(1:lx1,1:lx2,1:lx3)
      vs3(1:lx1,1:lx2,1:lx3,isp)=vs3(1:lx1,1:lx2,1:lx3,isp) &
                +muH(1:lx1,1:lx2,1:lx3,isp)/ns(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradpx(1:lx1,1:lx2,1:lx3) &
                -muP(1:lx1,1:lx2,1:lx3,isp)/ns(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradpy(1:lx1,1:lx2,1:lx3)    
    end do
  end if
end if

!> Gravitational drift terms (if required)
if (flaggravdrift) then
  if(flagswap/=1) then
    do isp=1,lsp
      vs2(1:lx1,1:lx2,1:lx3,isp)=vs2(1:lx1,1:lx2,1:lx3,isp)+ms(isp)/qs(isp)*(muP(:,:,:,isp)*g2-muH(:,:,:,isp)*g3)    !FIXME: +muH looks suspicious, I'm changing to (-)
      vs3(1:lx1,1:lx2,1:lx3,isp)=vs3(1:lx1,1:lx2,1:lx3,isp)+ms(isp)/qs(isp)*(muH(:,:,:,isp)*g2+muP(:,:,:,isp)*g3)
    end do
  else
    do isp=1,lsp
      vs2(1:lx1,1:lx2,1:lx3,isp)=vs2(1:lx1,1:lx2,1:lx3,isp)+ms(isp)/qs(isp)*(muP(:,:,:,isp)*g2+muH(:,:,:,isp)*g3)
      vs3(1:lx1,1:lx2,1:lx3,isp)=vs3(1:lx1,1:lx2,1:lx3,isp)+ms(isp)/qs(isp)*(-muH(:,:,:,isp)*g2+muP(:,:,:,isp)*g3)
    end do
  end if
end if


!! If it were appropriate this is how polarzations drifts could be computed.  However the particular quasistatic
!   model that we use explicitly omits this from the drift calculation which is then used in convective term in
!   polarization current
!    do isp=1,lsp
       !! To leading order the ion drifts do not include the polarization parts,
       !! otherwise it may mess up polarization convective term in the electrodynamics solver...
!      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE2Dt
!      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE3Dt
!    end do

end subroutine velocities


!>FIXME:  also needs the mobilities as input now due to pressure current terms.
subroutine potential_sourceterms(sigP,sigH,sigPgrav,sigHgrav,E02,E03,vn2,vn3,B1,muP,muH,ns,Ts,x, &
                                 flaggravdrift,flagdiamagnetic,srcterm)

!> Compute source terms (inhomogeneous terms) for the potential equation to be solved.  Both root and workers
!   should be able to use this routine

real(wp), dimension(:,:,:), intent(in) :: sigP,sigH,sigPgrav,sigHgrav
real(wp), dimension(:,:,:), intent(in) :: E02,E03,vn2,vn3
real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1     !ghost cells
real(wp), dimension(:,:,:,:), intent(in) :: muP,muH
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts     !ghost cells
type(curvmesh), intent(in) :: x
logical, intent(in) :: flaggravdrift
logical, intent(in) :: flagdiamagnetic
real(wp), dimension(:,:,:), intent(out) :: srcterm

real(wp), dimension(1:size(E02,1),1:size(E02,2),1:size(E02,3)) :: J1,J2,J3
real(wp), dimension(0:size(E02,1)+1,0:size(E02,2)+1,0:size(E02,3)+1) :: divtmp
!! one extra grid point on either end to facilitate derivatives
real(wp), dimension(-1:size(E02,1)+2,-1:size(E02,2)+2,-1:size(E02,3)+2) :: J1halo,J2halo,J3halo
!! haloing assumes existence of two ghost cells
integer :: lx1,lx2,lx3

real(wp), dimension(-1:size(E02,1)+2,-1:size(E02,2)+2,-1:size(E02,3)+2) :: pressure
real(wp), dimension(0:size(E02,1)+1,0:size(E02,2)+1,0:size(E02,3)+1) :: gradpx,gradpy


!> sizes from the conductivity coefficients
lx1=size(sigP,1)
lx2=size(sigP,2)
lx3=size(sigP,3)


!-------
!CONDUCTION CURRENT BACKGROUND SOURCE TERMS FOR POTENTIAL EQUATION. MUST COME AFTER CALL TO BC CODE.
J1 = 0                  
!! so this div is only perp components
J2 = 0 
J3 = 0     
!! for first current term zero everything out
call acc_perpconductioncurrents(sigP,sigH,E02,E03,J2,J3)     !background conduction currents only
if (debug .and. mpi_cfg%myid==0) print *, 'Workers has computed background field currents...'
call acc_perpwindcurrents(sigP,sigH,vn2,vn3,B1,J2,J3)
if (debug .and. mpi_cfg%myid==0) print *, 'Workers has computed wind currents...'
if (flagdiamagnetic) then
  call acc_pressurecurrents(muP,muH,ns,Ts,x,J2,J3)
  if (debug .and. mpi_cfg%myid==0) print *, 'Computed pressure currents...' 
end if
if (flaggravdrift) then
  call acc_perpgravcurrents(sigPgrav,sigHgrav,g2,g3,J2,J3)
  if (debug .and. mpi_cfg%myid==0) print *, 'Workers has computed gravitational currents...'
end if
!> FIXME:  need to add in pressure currents and then take divergence, also need mobilities for pressure currents.

J1halo(1:lx1,1:lx2,1:lx3)=J1     
!! temporary extended arrays to be populated with boundary data
J2halo(1:lx1,1:lx2,1:lx3)=J2
J3halo(1:lx1,1:lx2,1:lx3)=J3

call halo_pot(J1halo,tag%J1,x%flagper,.false.)
call halo_pot(J2halo,tag%J2,x%flagper,.false.)
call halo_pot(J3halo,tag%J3,x%flagper,.false.)

divtmp=div3D(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),J2halo(0:lx1+1,0:lx2+1,0:lx3+1), &
             J3halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
srcterm=divtmp(1:lx1,1:lx2,1:lx3)
!-------

end subroutine potential_sourceterms


subroutine acc_perpconductioncurrents(sigP,sigH,E2,E3,J2,J3)

!> ***Accumulate*** conduction currents into the variables J2,J3.  This
!    routine will not independently add background fields unless they are
!    already included in E2,3.  The currents are inout meaning that they
!    must be initialized to zero if you want only the conduction currents,
!    otherwise this routine just adds to whatever is already in J2,3.

real(wp), dimension(:,:,:), intent(in) :: sigP,sigH
real(wp), dimension(:,:,:), intent(in) :: E2,E3
real(wp), dimension(:,:,:), intent(inout) :: J2, J3


!> FIXME:  need to make consistent with velocity calculations flagswap /= 1
if (flagswap/=1) then
  J2=J2+sigP*E2-sigH*E3
  J3=J3+sigH*E2+sigP*E3
else
  J2=J2+sigP*E2+sigH*E3
  J3=J3-1*sigH*E2+sigP*E3
end if

end subroutine acc_perpconductioncurrents


subroutine acc_perpwindcurrents(sigP,sigH,vn2,vn3,B1,J2,J3)

!> ***Accumulate*** wind currents into the variables J2,J3.  See conduction currents
!    routine for additional caveats.

real(wp), dimension(:,:,:), intent(in) :: sigP,sigH
real(wp), dimension(:,:,:), intent(in) :: vn2,vn3
real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
real(wp), dimension(:,:,:), intent(inout) :: J2, J3

integer :: lx1,lx2,lx3
integer :: isp

!> sizes from the conductivities
lx1=size(sigP,1)
lx2=size(sigP,2)
lx3=size(sigP,3)

if (flagswap/=1) then     
  !! FIXME:  signs here require some explanation...  Perhaps add to formulation doc?
  J2=J2+sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
  J3=J3+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)-sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
else    !note signs on Pedersen terms change here (single cross product) but not Hall (double cross product generates cancelling sign changes)
  J2=J2-sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
  J3=J3+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)+sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
end if

end subroutine acc_perpwindcurrents


subroutine acc_pressurecurrents(muH,muP,ns,Ts,x,J2,J3)

!> ***Accumulate*** pressure currents into the variables J2,J3.  See conduction currents
!    routine for additional caveats.

real(wp), dimension(:,:,:,:), intent(in) :: muP,muH
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
type(curvmesh), intent(in) :: x
real(wp), dimension(:,:,:), intent(inout) :: J2, J3

real(wp), dimension(-1:size(J2,1)+2,-1:size(J2,2)+2,-1:size(J2,3)+2) :: pressure
real(wp), dimension(0:size(J2,1)+1,0:size(J2,2)+1,0:size(J2,3)+1) :: gradpx,gradpy

integer :: lx1,lx2,lx3,isp


!> sizes from the conductivities
lx1=size(J2,1)
lx2=size(J2,2)
lx3=size(J2,3)

if (flagswap/=1) then
  do isp=1,lsp
    pressure=ns(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)
    call halo_pot(pressure,tag%pressure,x%flagper,.false.)
    gradpx=grad3D2(pressure,x,0,lx1+1,0,lx2+1,0,lx3+1)
    gradpy=grad3D3(pressure,x,0,lx1+1,0,lx2+1,0,lx3+1)
    J2=J2+muP(:,:,:,isp)*gradpx(1:lx1,1:lx2,1:lx3)-muH(:,:,:,isp)*gradpy(1:lx1,1:lx2,1:lx3)
    J3=J3+muH(:,:,:,isp)*gradpx(1:lx1,1:lx2,1:lx3)+muP(:,:,:,isp)*gradpy(1:lx1,1:lx2,1:lx3)
  end do
else   
  do isp=1,lsp
    pressure=ns(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)
    call halo_pot(pressure,tag%pressure,x%flagper,.false.)
    gradpx=grad3D2(pressure,x,0,lx1+1,0,lx2+1,0,lx3+1)
    gradpy=grad3D3(pressure,x,0,lx1+1,0,lx2+1,0,lx3+1)
    J2=J2+muP(:,:,:,isp)*gradpx(1:lx1,1:lx2,1:lx3)+muH(:,:,:,isp)*gradpy(1:lx1,1:lx2,1:lx3)
    J3=J3-muH(:,:,:,isp)*gradpx(1:lx1,1:lx2,1:lx3)+muP(:,:,:,isp)*gradpy(1:lx1,1:lx2,1:lx3)
  end do
end if

end subroutine acc_pressurecurrents


subroutine acc_perpgravcurrents(sigPgrav,sigHgrav,g2,g3,J2,J3)

!> ***Accumulate*** gravitational currents into the variables J2,J3.  See conduction currents
!    routine for additional caveats.

real(wp), dimension(:,:,:), intent(in) :: sigPgrav,sigHgrav
real(wp), dimension(:,:,:), intent(in) :: g2,g3
real(wp), dimension(:,:,:), intent(inout) :: J2, J3


if (flagswap/=1) then
  J2=J2+sigPgrav*g2-sigHgrav*g3
  J3=J2+sigHgrav*g2+sigPgrav*g3
else
  J2=J2+sigPgrav*g2+sigHgrav*g3       !grav x2 current
  J3=J3-1*sigHgrav*g2+sigPgrav*g3     !grav x3 current
end if

end subroutine acc_perpgravcurrents


subroutine pot2perpfield(Phi,x,E2,E3)

!> computes electric field (perp components only) from a worker potential pattern.  Can
!   be called by either root or worker processes

real(wp), dimension(:,:,:), intent(in) :: Phi
type(curvmesh), intent(in) :: x
real(wp), dimension(:,:,:), intent(out) :: E2,E3

real(wp), dimension(0:size(Phi,1)+1,0:size(Phi,2)+1,0:size(Phi,3)+1) :: divtmp
!! one extra grid point on either end to facilitate derivatives
real(wp), dimension(-1:size(Phi,1)+2,-1:size(Phi,2)+2,-1:size(Phi,3)+2) :: J1halo,J2halo,J3halo
!! haloing assumes existence of two ghost cells
integer :: lx1,lx2,lx3


!> sizes from the mobility coefficients
lx1=size(Phi,1)
lx2=size(Phi,2)
lx3=size(Phi,3)

!CALCULATE PERP FIELDS FROM POTENTIAL
!      E20all=grad3D2(-1d0*Phi0all,dx2(1:lx2))
!! causes major memory leak. maybe from arithmetic statement argument?
!! Left here as a 'lesson learned' (or is it a gfortran bug...)
!      E30all=grad3D3(-1d0*Phi0all,dx3all(1:lx3all))

!COMPUTE THE 2 COMPONENT OF THE ELECTRIC FIELD
J1halo(1:lx1,1:lx2,1:lx3)=-1._wp*Phi
call halo_pot(J1halo,tag%J1,x%flagper,.true.)

divtmp=grad3D2(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
E2=divtmp(1:lx1,1:lx2,1:lx3)

!COMPUTE THE 3 COMPONENT OF THE ELECTRIC FIELD
J1halo(1:lx1,1:lx2,1:lx3)=-1._wp*Phi
call halo_pot(J1halo,tag%J1,x%flagper,.false.)

divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
E3=divtmp(1:lx1,1:lx2,1:lx3)
!--------

end subroutine pot2perpfield


subroutine get_BGEfields(x,E01,E02,E03)

type(curvmesh), intent(in) :: x
real(wp), dimension(:,:,:), intent(out) :: E01,E02,E03

!> routine to pull the background electric fields from the BCs modules and distribute
!   them to worker processes; called by both root and worker processes

real(wp), dimension(:,:,:), allocatable :: E01all,E02all,E03all    !space to pull the full dataset out of the module


if (mpi_cfg%myid==0) then
  allocate(E01all(lx1,lx2all,lx3all),E02all(lx1,lx2all,lx3all),E03all(lx1,lx2all,lx3all))
  E01all=0._wp
  call compute_rootBGEfields(x,E02all,E03all)

!  print*, 'E02all:  ',minval(E02all),maxval(E02all)
!  print*, 'E03all:  ',minval(E03all),maxval(E03all)

  call bcast_send(E01all,tag%E01,E01)
  call bcast_send(E02all,tag%E02,E02)
  call bcast_send(E03all,tag%E03,E03)
  deallocate(E01all,E02all,E03all)
else
  call bcast_recv(E01,tag%E01)
  call bcast_recv(E02,tag%E02)
  call bcast_recv(E03,tag%E03)
end if

!print*, mpi_cfg%myid,minval(E02),maxval(E02), shape(E02)
!print*, mpi_cfg%myid,minval(E03),maxval(E03), shape(E03)

end subroutine get_BGEfields


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

idleft=mpi_cfg%myid3-1
idright=mpi_cfg%myid3+1
iddown=mpi_cfg%myid2-1
idup=mpi_cfg%myid2+1

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
if (idup==mpi_cfg%lid2) then
  if (flagdegrade .and. lx2>1) then
    parmhalo(1:lx1,lx2+1,1:lx3)=2*parmhalo(1:lx1,lx2,1:lx3)-parmhalo(1:lx1,lx2-1,1:lx3)
  else
    parmhalo(1:lx1,lx2+1,1:lx3)=parmhalo(1:lx1,lx2,1:lx3)
  end if
end if
if (.not. flagper) then              
  !! musn't overwrite ghost cells if perioidc is chosen, only if aperiodic...
  if (idleft==-1) then
    parmhalo(1:lx1,1:lx2,0)=parmhalo(1:lx1,1:lx2,1)
  end if
  if (idright==mpi_cfg%lid3) then
    parmhalo(1:lx1,1:lx2,lx3+1)=parmhalo(1:lx1,1:lx2,lx3)
  end if
end if

end subroutine halo_pot

end module potential_comm

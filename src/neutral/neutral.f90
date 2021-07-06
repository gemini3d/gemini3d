module neutral

use phys_consts, only: wp, lnchem, pi, Re, debug
use grid, only: lx1, lx2, lx3
use meshobj, only: curvmesh
use timeutils, only : find_lastdate
use mpimod, only: mpi_cfg
use config, only: gemini_cfg

! also links MSIS from vendor/msis00/

implicit none (type, external)
private
public :: Tnmsis, neutral_atmos, make_dneu, clear_dneu, neutral_perturb, neutral_update, init_neutrals


interface ! atmos.f90
  module subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,v2grid,v3grid,nn,Tn,vn1,vn2,vn3, msis_version)
  integer, intent(in) :: ymd(3), msis_version
  real(wp), intent(in) :: UTsecd
  real(wp), dimension(:,:,:), intent(in) :: glat,glon,alt
  real(wp), intent(in) :: activ(3)
  real(wp), intent(in) :: v2grid,v3grid
  real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3),lnchem), intent(inout) :: nn
  !! intent(out)
  real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)), intent(inout) :: Tn
  !! intent(out)
  real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)), intent(inout) :: vn1,vn2,vn3
  !! intent(out)
  end subroutine neutral_atmos
end interface

interface ! perturb.f90
  module subroutine neutral_perturb(cfg,dt,dtneu,t,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
  type(gemini_cfg), intent(in) :: cfg
  real(wp), intent(in) :: dt,dtneu
  real(wp), intent(in) :: t
  integer, dimension(3), intent(in) :: ymd
  !! date for which we wish to calculate perturbations
  real(wp), intent(in) :: UTsec

  class(curvmesh), intent(inout) :: x
  !! grid structure  (inout becuase we want to be able to deallocate unit vectors once we are done with them)
  real(wp), intent(in) :: v2grid,v3grid
  real(wp), dimension(:,:,:,:), intent(inout) :: nn
  !! intent(out)
  !! neutral params interpolated to plasma grid at requested time
  real(wp), dimension(:,:,:), intent(inout) :: Tn,vn1,vn2,vn3
  !! intent(out)
  end subroutine neutral_perturb
end interface

!! ALL ARRAYS THAT FOLLOW ARE USED WHEN INCLUDING NEUTRAL PERTURBATIONS FROM ANOTHER MODEL
!! ARRAYS TO STORE THE NEUTRAL GRID INFORMATION
!! as long as the neutral module is in scope these persist and do not require a "save"; this variable only used by the axisymmetric interpolation
real(wp), dimension(:), allocatable, target, private :: rhon     !used for axisymmetric 2D simulations, aliased by pointers
real(wp), dimension(:), allocatable, target, private :: yn    !used in cartesian 2D and 3D interpolation
real(wp), dimension(:), allocatable, private :: zn
real(wp), dimension(:), allocatable, private :: xn    !for 3D cartesian interpolation
integer, private :: lrhon,lzn,lyn,lxn

!! STORAGE FOR NEUTRAL SIMULATION DATA.
! These will be singleton in the second dimension (longitude) in the case of 2D interpolation...
!! THESE ARE INCLUDED AS MODULE VARIATIONS TO AVOID HAVING TO REALLOCATE AND DEALLOCIATE EACH TIME WE NEED TO INTERP
real(wp), dimension(:,:,:), allocatable, private :: dnO,dnN2,dnO2,dvnrho,dvnz,dvnx,dTn

!!full grid parameters for root to store input from files.
real(wp), dimension(:), allocatable, private :: xnall
real(wp), dimension(:), allocatable, private :: ynall
integer, private :: lxnall,lynall

real(wp), dimension(:,:,:), allocatable, private :: dnOall,dnN2all,dnO2all,dvnrhoall,dvnzall,dvnxall,dTnall

!ARRAYS TO STORE NEUTRAL DATA THAT HAS BEEN INTERPOLATED
real(wp), dimension(:,:,:), allocatable, private :: dnOiprev,dnN2iprev,dnO2iprev,dvnrhoiprev,dvnziprev,dTniprev, &
                                                   dvn1iprev,dvn2iprev,dvn3iprev,dvnxiprev
real(wp), private :: tprev
integer, dimension(3), private :: ymdprev
!! time corresponding to "prev" interpolated data

real(wp), private :: UTsecprev
real(wp), dimension(:,:,:), allocatable, private :: dnOinext,dnN2inext,dnO2inext,dvnrhoinext,dvnzinext, &
                                                   dTninext,dvn1inext,dvn2inext,dvn3inext,dvnxinext
real(wp), private :: tnext
integer, dimension(3), private :: ymdnext
real(wp), private :: UTsecnext

!! data at current time level, (centered in time between current time step and next)
real(wp), dimension(:,:,:), allocatable, protected :: dnOinow,dnN2inow,dnO2inow,dTninow,dvn1inow,dvn2inow,dvn3inow

!SPACE TO STORE PROJECTION FACTORS (rotate from magnetic UEN to curv. dipole...)
real(wp), dimension(:,:,:), allocatable, private :: proj_erhop_e1,proj_ezp_e1,proj_erhop_e2,proj_ezp_e2,proj_erhop_e3,proj_ezp_e3    !these projections are used in the axisymmetric interpolation
real(wp), dimension(:,:,:), allocatable, private :: proj_eyp_e1,proj_eyp_e2,proj_eyp_e3    !these are for Cartesian projections
real(wp), dimension(:,:,:), allocatable, private :: proj_exp_e1,proj_exp_e2,proj_exp_e3

!PLASMA GRID ZI AND RHOI LOCATIONS FOR INTERPOLATIONS
real(wp), dimension(:), allocatable, private :: zi,xi    !this is to be a flat listing of sites on the, rhoi only used in axisymmetric and yi only in cartesian
real(wp), dimension(:), allocatable, target, private :: yi,rhoi

!USED FOR 3D INTERPOLATION WHERE WORKER DIVISIONS ARE COMPLICATED (note that the first dim starts at zero so it matches mpi ID)
real(wp), dimension(:,:), private, allocatable :: extents    !roots array that is used to store min/max x,y,z of each works
integer, dimension(:,:), private, allocatable :: indx       !roots array that contain indices for each workers needed piece of the neutral data
integer, dimension(:,:), private, allocatable :: slabsizes

!! BASE MSIS ATMOSPHERIC STATE ON WHICH TO APPLY PERTURBATIONS
real(wp), dimension(:,:,:,:), allocatable, protected :: nnmsis
real(wp), dimension(:,:,:), allocatable, protected :: Tnmsis
real(wp), dimension(:,:,:), allocatable, protected :: vn1base,vn2base,vn3base

!! projection factors for converting vectors mag->geo; e.g. defining rotation matrix from geographic coords into
real(wp), dimension(:,:,:), allocatable :: proj_ealt_e1,proj_ealt_e2,proj_ealt_e3
real(wp), dimension(:,:,:), allocatable :: proj_eglat_e1,proj_eglat_e2,proj_eglat_e3
real(wp), dimension(:,:,:), allocatable :: proj_eglon_e1,proj_eglon_e2,proj_eglon_e3

contains


subroutine init_neutrals(dt,t,cfg,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)

!> initializes neutral atmosphere by:
!    1)  allocating storage space
!    2)  establishing initial background
!    3)  priming file input so that we have an initial perturbed state to start from (necessary for restart)

real(wp), intent(in) :: dt,t
type(gemini_cfg), intent(in) :: cfg
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
class(curvmesh), intent(inout) :: x
!! unit vecs may be deallocated after first setup
real(wp), intent(in) :: v2grid,v3grid
real(wp), dimension(:,:,:,:), intent(inout) :: nn
!! intent(out)
real(wp), dimension(:,:,:), intent(inout) :: Tn
!! intent(out)
real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
!! intent(out)
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp

real(wp) :: tstart,tfin


!! allocation neutral module scope variables so there is space to store all the file input and do interpolations
call make_dneu()

!! call msis to get an initial neutral background atmosphere
if (mpi_cfg%myid == 0) call cpu_time(tstart)
call neutral_atmos(ymd,UTsec,x%glat,x%glon,x%alt,cfg%activ,v2grid,v3grid,nn,Tn,vn1,vn2,vn3, cfg%msis_version)
if (mpi_cfg%myid == 0) then
  call cpu_time(tfin)
  print *, 'Initial neutral background at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
end if

if (cfg%flagdneu==1) then
  !! find the last input data preceding the milestone/initial condition that we start with
  call find_lastdate(cfg%ymd0,cfg%UTsec0,ymd,UTsec,cfg%dtneu,ymdtmp,UTsectmp)

  !! Loads the neutral input file corresponding to the "first" time step of the simulation to prevent the first interpolant
  !  from being zero and causing issues with restart simulations.  I.e. make sure the neutral buffers are primed for restart
  !  This requires us to load file input twice, once corresponding to the initial frame and once for the "first, next" frame.
  tprev=UTsectmp-UTsec-2*cfg%dtneu
  tnext=tprev+cfg%dtneu
  if (mpi_cfg%myid==0) print*, '!!!Attempting initial load of neutral dynamics files!!!' // &
                           ' This is a workaround to insure compatibility with restarts...',ymdtmp,UTsectmp
  !! We essentially are loading up the data corresponding to halfway betwween -dtneu and t0 (zero).  This will load
  !   two time levels back so when tprev is incremented twice it will be the true tprev corresponding to first time step
  call neutral_perturb(cfg,dt,cfg%dtneu,tnext+cfg%dtneu/2,ymdtmp,UTsectmp-cfg%dtneu, &
                        x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)  !abs time arg to be < 0

  if (mpi_cfg%myid==0) print*, 'Now loading initial next file for neutral perturbations...'
  !! Now compute perturbations for the present time (zero), this moves the primed variables in next into prev and then
  !  loads up a current state so that we get a proper interpolation for the first time step.
  call neutral_perturb(cfg,dt,cfg%dtneu,0._wp,ymdtmp,UTsectmp,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)    !t-dt so we land exactly on start time
end if
end subroutine init_neutrals


subroutine neutral_update(nn,Tn,vn1,vn2,vn3)

!! adds stored base and perturbation neutral atmospheric parameters
!!  these are module-scope parameters so not needed as input

real(wp), dimension(:,:,:,:), intent(inout) :: nn
!! intent(out)
real(wp), dimension(:,:,:), intent(inout) :: Tn
!! intent(out)
real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
!! intent(out)


!> background neutral parameters
nn=nnmsis
Tn=Tnmsis
vn1=vn1base
vn2=vn2base
vn3=vn3base

!> perturbations, if used
if (allocated(zn)) then
  nn(:,:,:,1)=nn(:,:,:,1)+dnOinow
  nn(:,:,:,2)=nn(:,:,:,2)+dnN2inow
  nn(:,:,:,3)=nn(:,:,:,3)+dnO2inow
  nn(:,:,:,1)=max(nn(:,:,:,1),1._wp)
  nn(:,:,:,2)=max(nn(:,:,:,2),1._wp)
  nn(:,:,:,3)=max(nn(:,:,:,3),1._wp)
  !! note we are not adjusting derived densities like NO since it's not clear how they may be related to major
  !! species perturbations.

  Tn=Tn+dTninow
  Tn=max(Tn,50._wp)

  vn1=vn1+dvn1inow
  vn2=vn2+dvn2inow
  vn3=vn3+dvn3inow
end if

end subroutine neutral_update


!> rotate winds from geographic to model native coordinate system (x1,x2,x3)
subroutine rotate_geo2native(vnalt,vnglat,vnglon,x,vn1,vn2,vn3)
  real(wp), dimension(:,:,:), intent(in) :: vnalt,vnglat,vnglon
  class(curvmesh), intent(in) :: x
  real(wp), dimension(:,:,:), intent(out) :: vn1,vn2,vn3
  real(wp), dimension(1:size(vnalt,1),1:size(vnalt,2),1:size(vnalt,3)) :: ealt,eglat,eglon
  integer :: lx1,lx2,lx3

  !> if first time called then allocate space for projections and compute
  if (.not. allocated(proj_ealt_e1) then
    call x%calc_unitvec_geo(ealt,eglat,eglon)

    lx1=size(vnalt,1); lx2=size(vnalt,2); lx3=size(vnalt,3);
    allocate(proj_ealt_e1(lx1,lx2,lx3),proj_eglat_e1(lx1,lx2,lx3),proj_eglon_e1(lx1,lx2,lx3))
    allocate(proj_ealt_e2(lx1,lx2,lx3),proj_eglat_e2(lx1,lx2,lx3),proj_eglon_e2(lx1,lx2,lx3))
    allocate(proj_ealt_e3(lx1,lx2,lx3),proj_eglat_e3(lx1,lx2,lx3),proj_eglon_e3(lx1,lx2,lx3))

    proj_ealt_e1=sum(ealt*x%e1,4)
    proj_eglat_e1=sum(eglat*x%e1,4)
    proj_eglon_e1=sum(eglon*x%e1,4)
    proj_ealt_e2=sum(ealt*x%e2,4)
    proj_eglat_e2=sum(eglat*x%e2,4)
    proj_eglon_e2=sum(eglon*x%e2,4)
    proj_ealt_e3=sum(ealt*x%e3,4)
    proj_eglat_e3=sum(eglat*x%e3,4)
    proj_eglon_e3=sum(eglon*x%e3,4)
  end if

  !> rotate vectors into model native coordinate system
  vn1=vnalt*proj_ealt_e1+vnglat*proj_eglat_e1+vnglon*proj_eglon_e1
  vn2=vnalt*proj_ealt_e2+vnglat*proj_eglat_e2+vnglon*proj_eglon_e2
  vn3=vnalt*proj_ealt_e3+vnglat*proj_eglat_e3+vnglon*proj_eglon_e3
end subroutine rotate_geo2native


subroutine make_dneu()
!ZZZ - could make this take in type of neutral interpolation and do allocations accordingly

!allocate and compute plasma grid z,rho locations and space to save neutral perturbation variables and projection factors
allocate(zi(lx1*lx2*lx3),rhoi(lx1*lx2*lx3))
allocate(yi(lx1*lx2*lx3))
allocate(xi(lx1*lx2*lx3))
allocate(proj_erhop_e1(lx1,lx2,lx3),proj_ezp_e1(lx1,lx2,lx3),proj_erhop_e2(lx1,lx2,lx3),proj_ezp_e2(lx1,lx2,lx3), &
         proj_erhop_e3(lx1,lx2,lx3),proj_ezp_e3(lx1,lx2,lx3))
allocate(proj_eyp_e1(lx1,lx2,lx3),proj_eyp_e2(lx1,lx2,lx3),proj_eyp_e3(lx1,lx2,lx3))
allocate(proj_exp_e1(lx1,lx2,lx3),proj_exp_e2(lx1,lx2,lx3),proj_exp_e3(lx1,lx2,lx3))
allocate(dnOiprev(lx1,lx2,lx3),dnN2iprev(lx1,lx2,lx3),dnO2iprev(lx1,lx2,lx3),dvnrhoiprev(lx1,lx2,lx3), &
         dvnziprev(lx1,lx2,lx3),dTniprev(lx1,lx2,lx3),dvn1iprev(lx1,lx2,lx3),dvn2iprev(lx1,lx2,lx3), &
         dvn3iprev(lx1,lx2,lx3))
allocate(dvnxiprev(lx1,lx2,lx3))
allocate(dnOinext(lx1,lx2,lx3),dnN2inext(lx1,lx2,lx3),dnO2inext(lx1,lx2,lx3),dvnrhoinext(lx1,lx2,lx3), &
         dvnzinext(lx1,lx2,lx3),dTninext(lx1,lx2,lx3),dvn1inext(lx1,lx2,lx3),dvn2inext(lx1,lx2,lx3), &
         dvn3inext(lx1,lx2,lx3))
allocate(dvnxinext(lx1,lx2,lx3))
allocate(nnmsis(lx1,lx2,lx3,lnchem),Tnmsis(lx1,lx2,lx3),vn1base(lx1,lx2,lx3),vn2base(lx1,lx2,lx3),vn3base(lx1,lx2,lx3))
allocate(dnOinow(lx1,lx2,lx3),dnN2inow(lx1,lx2,lx3),dnO2inow(lx1,lx2,lx3),dvn1inow(lx1,lx2,lx3),dvn2inow(lx1,lx2,lx3), &
           dvn3inow(lx1,lx2,lx3), dTninow(lx1,lx2,lx3))

!start everyone out at zero
zi = 0
rhoi = 0
yi = 0
xi = 0
proj_erhop_e1 = 0
proj_ezp_e1 = 0
proj_erhop_e2 = 0
proj_ezp_e2 = 0
proj_erhop_e3 = 0
proj_ezp_e3 = 0
proj_eyp_e1 = 0
proj_eyp_e2 = 0
proj_eyp_e3 = 0
proj_exp_e1 = 0
proj_exp_e2 = 0
proj_exp_e3 = 0
dnOiprev = 0
dnN2iprev = 0
dnO2iprev = 0
dTniprev = 0
dvnrhoiprev = 0
dvnziprev = 0
dvn1iprev = 0
dvn2iprev = 0
dvn3iprev = 0
dvnxiprev = 0
dnOinext = 0
dnN2inext = 0
dnO2inext = 0
dTninext = 0
dvnrhoinext = 0
dvnzinext = 0
dvn1inext = 0
dvn2inext = 0
dvn3inext = 0
dvnxinext = 0
nnmsis = 0
Tnmsis = 0
vn1base = 0
vn2base = 0
vn3base = 0
dnOinow = 0
dnN2inow = 0
dnO2inow = 0
dTninow = 0
dvn1inow = 0
dvn2inow = 0
dvn3inow = 0

!now initialize some module variables
tprev = 0
tnext = 0

end subroutine make_dneu


subroutine clear_dneu

!stuff allocated at beginning of program
deallocate(zi,rhoi)
deallocate(yi)
deallocate(proj_erhop_e1,proj_ezp_e1,proj_erhop_e2,proj_ezp_e2, &
         proj_erhop_e3,proj_ezp_e3)
deallocate(proj_eyp_e1,proj_eyp_e2,proj_eyp_e3)
deallocate(dnOiprev,dnN2iprev,dnO2iprev,dvnrhoiprev,dvnziprev,dTniprev,dvn1iprev,dvn2iprev,dvn3iprev)
deallocate(dnOinext,dnN2inext,dnO2inext,dvnrhoinext,dvnzinext,dTninext,dvn1inext,dvn2inext,dvn3inext)
deallocate(nnmsis,Tnmsis,vn1base,vn2base,vn3base)
deallocate(dnOinow,dnN2inow,dnO2inow,dTninow,dvn1inow,dvn2inow,dvn3inow)

!check whether any other module variables were allocated and deallocate accordingly
if (allocated(zn) ) then    !if one is allocated, then they all are
  deallocate(zn)
  deallocate(dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
end if
if (allocated(rhon)) then
  deallocate(rhon)
end if
if (allocated(yn)) then
  deallocate(yn)
end if
if (allocated(extents)) then
  deallocate(extents,indx,slabsizes)
end if
if (allocated(dvnx)) then
  deallocate(dvnx)
end if
if (allocated(xn)) then
  deallocate(xn)
end if
if (allocated(xnall)) then
  deallocate(xnall,ynall)
  deallocate(dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
end if
if (allocated(proj_ealt_e1) then
  deallocate(proj_ealt_e1,proj_eglat_e1,proj_eglon_e1)
  deallocate(proj_ealt_e2,proj_eglat_e2,proj_eglon_e2)
  deallocate(proj_ealt_e3,proj_eglat_e3,proj_eglon_e3)
end if

end subroutine clear_dneu

end module neutral

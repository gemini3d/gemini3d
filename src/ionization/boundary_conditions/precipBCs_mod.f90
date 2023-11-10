module precipBCs_mod

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

use reader, only: get_simsize2, get_precip, get_grid2
use phys_consts, only: pi,wp, debug
use grid, only : lx1,lx2,lx3
use meshobj, only: curvmesh
use interpolation, only : interp1,interp2
use timeutils, only : dateinc, date_filename, find_lastdate
use mpimod, only: mpi_realprec, mpi_cfg, tag=>gemini_mpi
use gemini3d_config, only: gemini_cfg
use precipdataobj, only: precipdata

use mpi_f08, only : mpi_integer, mpi_comm_world, mpi_status_ignore

implicit none (type, external)
private
public :: precipBCs_fileinput, precipBCs, init_precipinput

contains
  !> initialize variables to hold input file precipitation information, must be called by all workers at the same time
  subroutine init_precipinput(dt,cfg,ymd,UTsec,x,eprecip)
    real(wp), intent(in) :: dt
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    type(precipdata), intent(inout) :: eprecip

    if (cfg%flagprecfile==1) then    !all workers must have this info
      call eprecip%init(cfg,cfg%precdir,x,dt,cfg%dtprec,ymd,UTsec)
    end if
  end subroutine init_precipinput


  !> get latest file input precipitation information
  subroutine precipBCs_fileinput(dtmodel,t,cfg,ymd,UTsec,x,W0,PhiWmWm2,eprecip)
    real(wp), intent(in) :: dtmodel
    real(wp), intent(in) :: t
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    !! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:), intent(inout) :: W0,PhiWmWm2
    !! intent(out)
    !! last dimension is the number of particle populations
    class(curvmesh), intent(in) :: x
    type(precipdata), intent(inout) :: eprecip
    integer :: ix2,ix3

    ! background precipitation from config.nml file
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        W0(ix2,ix3,1)=cfg%W0BG
        PhiWmWm2(ix2,ix3,1)=cfg%PhiWBG
      end do
    end do

    ! disturbance precipitation from file input
    call eprecip%update(cfg,dtmodel,t,x,ymd,UTsec)

    ! set output arrays; note that this is making a copy of the data stored the precipdata object.  We
    !   will assume for now that this doesn't incur too much memory overhead
    W0(:,:,2)=eprecip%E0pinow(:,:)
    PhiWmWm2(:,:,2)=eprecip%Qpinow(:,:)

    !print*, minval(W0(:,:,2)),maxval(W0(:,:,2))
    !print*, minval(PhiWmWm2(:,:,2)),maxval(PhiWmWm2(:,:,2))
  end subroutine precipBCs_fileinput


  !> This is the default subroutine that is called for electron precipitation if file input is not used.
  subroutine precipBCs(cfg,W0,PhiWmWm2)
    !------------------------------------------------------------
    !-------LOAD UP ARRAYS CONTAINING TOP BOUNDARY CHAR. ENERGY
    !-------AND TOTAL ENERGY FLUX.  GRID VARIABLES INCLUDE
    !-------GHOST CELLS
    !------------------------------------------------------------

    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:), intent(inout) :: W0,PhiWmWm2
    !! intent(out)
    real(wp) :: W0pk,PhiWpk
    real(wp) :: sigx2,sigx3,x30amp,varc,meanx2,sigt,meant
    integer :: ix2,ix3,lx2,lx3,lprec


    lx2=size(W0,1)
    lx3=size(W0,2)
    lprec=size(W0,3)    !assumed to be 2 in this subroutine


    !BACKGROUND PRECIPITATION
    W0pk = cfg%W0BG
    PhiWpk=cfg%PhiWBG
    !PhiWpk = 1e-3_wp
    do ix3=1,lx3
      do ix2=1,lx2
        W0(ix2,ix3,1)=W0pk
        PhiWmWm2(ix2,ix3,1)=PhiWpk
      end do
    end do


    !PARAMETERS FOR DISTURBANCE PRECIPITATION
    W0pk = 100
    !      sigW0x3=100e3_wp
    !      meanW0x3=0
    PhiWpk = 0
    !      PhiWpk=1e-5_wp    !successful grad-drift attempts
    !      PhiWpk=1e-4_wp   !Swoboda blur testing
    !      PhiWpk=0.05_wp    !testing of convergent Hall drifts
    !      PhiWpk=5._wp
    !      sigPhiWx3=100e3_wp
    !      meanPhiWx3=0

    !      W0pk=0.3e3_wp
    !      sigW0x3=100e3_wp
    !      meanW0x3=0
    !      PhiWpk=2._wp
    !      sigPhiWx3=100e3_wp
    !      meanPhiWx3=0

    sigx2 = 50e3_wp
    meanx2 = 0
    !    sigx3=10e3_wp
    sigx3 = 25e3_wp
    meant = 900
    sigt = 450
    x30amp= 0
    varc = 200

    !DISTURBANCE ELECTRON PRECIPITATION PATTERN
    do ix3=1,lx3
      do ix2=1,lx2
        W0(ix2,ix3,2) = W0pk
        PhiWmWm2(ix2,ix3,2) = PhiWpk
      end do
    end do
  end subroutine precipBCs
end module precipBCs_mod

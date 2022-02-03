!> mpi-related gemini functionality; all procedures must be bind(C) to allow C/C++ main programs
module gemini3d_mpi

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use phys_consts, only: wp,debug
use mpimod, only: mpi_manualgrid, process_grid_auto, mpi_cfg, mpibreakdown
use meshobj, only: curvmesh
use config, only: gemini_cfg
use io, only: output_plasma,output_aur,find_milestone,input_plasma,create_outdir,create_outdir_aur
use potential_comm, only: get_BGEfields,velocities
use grid_mpi, only: grid_drift
use collisions, only: conductivities
use gemini3d, only: cfg,x

implicit none (type, external)
private
public :: init_procgrid, outdir_fullgridvaralloc, get_initial_state, BGfield_Lagrangian, check_dryrun, check_fileoutput,  &
            get_initial_drifts

contains
  !> establish gemini process grid
  subroutine init_procgrid(lx2all,lx3all,lid2in,lid3in) 
    integer, intent(in) :: lx2all,lx3all,lid2in,lid3in
    if (lid2in==-1) then
      call process_grid_auto(lx2all, lx3all)
      !! grid_size defines lx2all and lx3all
    else
      call mpi_manualgrid(lx2all, lx3all, lid2in, lid3in)
    endif
    print '(A, I0, A1, I0)', 'process grid (Number MPI processes) x2, x3:  ',mpi_cfg%lid2, ' ', mpi_cfg%lid3
    print '(A, I0, A, I0, A1, I0)', 'Process:',mpi_cfg%myid,' at process grid location: ',mpi_cfg%myid2,' ',mpi_cfg%myid3
  end subroutine init_procgrid

    !> add in background field, accounting for whether the user specified a lagrangian grid
!  module subroutine BGfield_Lagrangian(cfg,x,v2grid,v3grid,E1,E2,E3)
!    type(gemini_cfg), intent(in) :: cfg
!    class(curvmesh), intent(in) :: x
!    real(wp), intent(inout) :: v2grid,v3grid
!    real(wp), dimension(:,:,:), intent(inout) :: E1,E2,E3
!  end subroutine BGfield_Lagrangian
  subroutine BGfield_Lagrangian(v2grid,v3grid,E1,E2,E3) 
    real(wp), intent(inout) :: v2grid,v3grid
    real(wp), dimension(:,:,:), intent(inout) :: E1,E2,E3
    real(wp), dimension(:,:,:), allocatable :: E01,E02,E03
    integer :: lx1,lx2,lx3

    lx1=size(E2,1); lx2=size(E2,2); lx3=size(E3,3);
    allocate(E01(lx1,lx2,lx3),E02(lx1,lx2,lx3),E03(lx1,lx2,lx3))
    E01=0; E02=0; E03=0;
    if (cfg%flagE0file==1) then
      call get_BGEfields(x,E01,E02,E03)
    end if
    if (cfg%flaglagrangian) then    ! Lagrangian (moving) grid; compute from input background electric fields
      call grid_drift(x,E02,E03,v2grid,v3grid)
      if (mpi_cfg%myid==0) print*, mpi_cfg%myid,' using Lagrangian grid moving at:  ',v2grid,v3grid
    else                            ! stationary grid
      v2grid = 0
      v3grid = 0
      E1 = E1 + E01    ! FIXME: this is before dist fields are computed???
      E2 = E2 + E02
      E3 = E3 + E03
    end if
    deallocate(E01,E02,E03)
  end subroutine BGfield_Lagrangian

  !> Compute initial perp drifts
!  module subroutine get_initial_drifts(cfg,x,nn,Tn,vn1,vn2,vn3,ns,Ts,vs1,vs2,vs3,B1,E2,E3)
!    type(gemini_cfg), intent(in) :: cfg
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:,:), intent(in) :: nn
!    real(wp), dimension(:,:,:), intent(in) :: Tn,vn1,vn2,vn3
!    real(wp), dimension(:,:,:,:), intent(in) :: ns,Ts,vs1
!    real(wp), dimension(:,:,:,:), intent(inout) :: vs2,vs3
!    real(wp), dimension(:,:,:), intent(in) :: B1
!    real(wp), dimension(:,:,:), intent(in) :: E2,E3
!  end subroutine get_initial_drifts
  subroutine get_initial_drifts(nn,Tn,vn1,vn2,vn3,ns,Ts,vs1,vs2,vs3,B1,E2,E3) 
    real(wp), dimension(:,:,:,:), intent(in) :: nn
    real(wp), dimension(:,:,:), intent(in) :: Tn,vn1,vn2,vn3
    real(wp), dimension(:,:,:,:), intent(in) :: ns,Ts,vs1
    real(wp), dimension(:,:,:,:), intent(inout) :: vs2,vs3
    real(wp), dimension(:,:,:), intent(in) :: B1
    real(wp), dimension(:,:,:), intent(in) :: E2,E3
    real(wp), dimension(:,:,:), allocatable :: sig0,sigP,sigH,sigPgrav,sigHgrav
    real(wp), dimension(:,:,:,:), allocatable :: muP,muH,nusn
    integer :: lx1,lx2,lx3,lsp
    
    lx1=x%lx1; lx2=x%lx2; lx3=x%lx3; lsp=size(ns,4);
    allocate(sig0(lx1,lx2,lx3),sigP(lx1,lx2,lx3),sigH(lx1,lx2,lx3),sigPgrav(lx1,lx2,lx3),sigHgrav(lx1,lx2,lx3))
    allocate(muP(lx1,lx2,lx3,lsp),muH(lx1,lx2,lx3,lsp),nusn(lx1,lx2,lx3,lsp))
    call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
    call velocities(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)
    deallocate(sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
  end subroutine get_initial_drifts

  !> Create output directories and allocate root-only variables
  subroutine outdir_fullgridvaralloc(Phiall,lx1,lx2all,lx3all)
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall
    integer, intent(in) :: lx1,lx2all,lx3all
    !> create a place, if necessary, for output datafiles 
    if (mpi_cfg%myid==0) then
      call create_outdir(cfg)
      if (cfg%flagglow /= 0) call create_outdir_aur(cfg%outdir)
    end if  

    !> fullgrid variable allocations only needed for the potential variable
    if (mpi_cfg%myid==0) then
      allocate(Phiall(lx1,lx2all,lx3all))
    end if
  end subroutine outdir_fullgridvaralloc


  !> Determine whether we are restarting vs. starting from a user-specified state.  Note that this uses mpi right now but
  !    with Michael's h5fortran-mpi library this call will be executed by all workers
  subroutine get_initial_state(ns,vs1,Ts,Phi,Phiall,UTsec,ymd,tdur)
    real(wp), dimension(:,:,:,:), intent(inout) :: ns,vs1,Ts
    real(wp), dimension(:,:,:), intent(inout) :: Phi,Phiall
    real(wp), intent(inout) :: UTsec
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: tdur
  !end subroutine get_initial_state
    integer, dimension(3) :: ymdtmp
    real(wp) :: UTsectmp,ttmp
    character(:), allocatable :: filetmp

    call find_milestone(cfg, ttmp, ymdtmp, UTsectmp, filetmp)
    if ( ttmp > 0 ) then
      !! restart scenario
      if (mpi_cfg%myid==0) then
        print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print*, '! Restarting simulation from time:  ',ymdtmp,UTsectmp
        print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      end if
    
      !! Set start variables accordingly and read in the milestone
      UTsec=UTsectmp
      ymd=ymdtmp
      tdur=cfg%tdur-ttmp    ! subtract off time that has elapsed to milestone
      if (mpi_cfg%myid==0) then
        print*, 'Treating the following file as initial conditions:  ',filetmp
        print*, ' full duration:  ',cfg%tdur,'; remaining simulation time:  ',tdur
      end if
    
      if (tdur <= 1e-6_wp .and. mpi_cfg%myid==0) error stop 'Cannot restart simulation from the final time step!'
    
      cfg%tdur=tdur         ! just to insure consistency
      call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,filetmp,ns,vs1,Ts,Phi,Phiall)
    else !! start at the beginning
      UTsec = cfg%UTsec0
      ymd = cfg%ymd0
      tdur = cfg%tdur
    
      if (tdur <= 1e-6_wp .and. mpi_cfg%myid==0) error stop 'Simulation is of zero time duration'
      call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,cfg%indatfile,ns,vs1,Ts,Phi,Phiall)
    end if
  end subroutine get_initial_state


  !> see if we need to perform an output
  subroutine check_fileoutput(t,tout,tglowout,tmilestone,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3,iver)
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: tout,tglowout,tmilestone
    integer, intent(inout) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall
    real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3
    real(wp), dimension(:,:,:), intent(in) :: iver
!  end subroutine
    real(wp) :: tstart,tfin


    if (abs(t-tout) < 1d-5) then
      tout = tout + cfg%dtout
      if (cfg%nooutput ) then
        if (mpi_cfg%myid==0) write(stderr,*) 'WARNING: skipping file output at sim time (sec)',t
        return
      endif
      !! close enough to warrant an output now...
      if (mpi_cfg%myid==0 .and. debug) call cpu_time(tstart)
  
      !! We may need to adjust flagoutput if we are hitting a milestone
      flagoutput=cfg%flagoutput
      if (cfg%mcadence>0 .and. abs(t-tmilestone) < 1d-5) then
        flagoutput=1    !force a full output at the milestone
        call output_plasma(cfg%outdir,flagoutput,ymd, &
          UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, &
          cfg%out_format)
        tmilestone = t + cfg%dtout * cfg%mcadence
        if(mpi_cfg%myid==0) print*, 'Milestone output triggered.'
      else
        call output_plasma(cfg%outdir,flagoutput,ymd, &
          UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, &
          cfg%out_format)
      end if
      if (mpi_cfg%myid==0 .and. debug) then
        call cpu_time(tfin)
        print *, 'Plasma output done for time step:  ',t,' in cpu_time of:  ',tfin-tstart
      endif
    end if
  
    !> GLOW file output
    if ((cfg%flagglow /= 0) .and. (abs(t-tglowout) < 1d-5)) then !same as plasma output
      call cpu_time(tstart)
      call output_aur(cfg%outdir, cfg%flagglow, ymd, UTsec, iver, cfg%out_format)
      if (mpi_cfg%myid==0) then
        call cpu_time(tfin)
        print *, 'Auroral output done for time step:  ',t,' in cpu_time of: ',tfin-tstart
      end if
      tglowout = tglowout + cfg%dtglowout
    end if
  end subroutine check_fileoutput


  !> check whether user called for a dryrun and end the program if so
  subroutine check_dryrun
    character(8) :: date
    character(10) :: time
    integer :: ierr

    if (cfg%dryrun) then
      ierr = mpibreakdown()
      if (ierr /= 0) error stop 'Gemini dry run MPI shutdown failure'
      call date_and_time(date,time)
      print '(/,A)', 'DONE: ' // date(1:4) // '-' // date(5:6) // '-' // date(7:8) // 'T' &
        // time(1:2) // ':' // time(3:4) // ':' // time(5:)
      stop "OK: Gemini dry run"
    endif
  end subroutine check_dryrun 
end module gemini3d_mpi

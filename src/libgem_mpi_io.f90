!> io related subroutines for the gemini3d_mpi library
submodule (gemini3d_mpi) libgem_mpi_io

implicit none (type, external)

contains
  !> Create output directories and allocate root-only variables
  !module subroutine outdir_fullgridvaralloc(cfg,Phiall,lx1,lx2all,lx3all)
  !  type(gemini_cfg), intent(in) :: cfg
  !  real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall
  !  integer, intent(in) :: lx1,lx2all,lx3all
  !end subroutine outdir_fullgridvaralloc
  module procedure outdir_fullgridvaralloc
    !> create a place, if necessary, for output datafiles
    if (mpi_cfg%myid==0) then
      call create_outdir(cfg)
    end if

    !> fullgrid variable allocations only needed for the potential variable
    if (mpi_cfg%myid==0) then
      allocate(Phiall(lx1,lx2all,lx3all))
    end if
  end procedure outdir_fullgridvaralloc


  !> Determine whether we are restarting vs. starting from a user-specified state.  Note that this uses mpi right now but
  !    with Michael's h5fortran-mpi library this call will be executed by all workers
  !module subroutine get_initial_state(cfg,x,ns,vs1,Ts,Phi,Phiall,UTsec,ymd,tdur)
  !  type(gemini_cfg), intent(inout) :: cfg
  !  class(curvmesh), intent(in) :: x
  !  real(wp), dimension(:,:,:,:), intent(inout) :: ns,vs1,Ts
  !  real(wp), dimension(:,:,:), intent(inout) :: Phi,Phiall
  !  real(wp), intent(inout) :: UTsec
  !  integer, dimension(3), intent(inout) :: ymd
  !  real(wp), intent(inout) :: tdur
  !end subroutine get_initial_state
  module procedure get_initial_state
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
  end procedure get_initial_state


  !> see if we need to perform an output
!  module subroutine check_fileoutput(t,tout,tglowout,tmilestone,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3,iver)
!    real(wp), intent(in) :: t
!    real(wp), intent(inout) :: tout,tglowout,tmilestone
!    type(gemini_cfg), intent(in) :: cfg
!    integer, intent(inout) :: flagoutput
!    integer, dimension(3), intent(in) :: ymd
!    real(wp), intent(in) :: UTsec
!    real(wp), dimension(:,:,:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
!    real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall
!    real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3
!    real(wp), dimension(:,:,:), intent(in) :: iver
!  end subroutine
  module procedure check_fileoutput
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
  end procedure check_fileoutput


  !> check whether user called for a dryrun and end the program if so
  module procedure check_dryrun
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
  end procedure check_dryrun
end submodule libgem_mpi_io

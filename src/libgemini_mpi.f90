!> mpi-related gemini functionality
module gemini3d_mpi

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: iso_c_binding, only : c_ptr
use phys_consts, only: wp,debug
use mpimod, only: mpi_manualgrid, process_grid_auto, mpi_cfg, mpibreakdown, mpisetup
use meshobj, only: curvmesh
use config, only: gemini_cfg
use io, only: output_plasma,output_aur,find_milestone,input_plasma,create_outdir
use potential_comm, only: get_BGEfields,velocities
use grid_mpi, only: grid_drift, read_grid
use collisions, only: conductivities
use potentialBCs_mumps, only: init_Efieldinput
use potential_comm,only : pot2perpfield
use neutral_perturbations, only: init_neutralperturb,neutral_denstemp_update,neutral_wind_update,neutral_perturb
use temporal, only : dt_comm
use sanity_check, only : check_finite_pertub, check_finite_output
use potential_comm,only : electrodynamics
use advec_mpi, only: set_global_boundaries_allspec, halo_interface_vels_allspec
use multifluid_mpi, only: halo_allparams
use sources_mpi, only: RK2_prep_mpi_allspec
use ionization_mpi, only: get_gavg_Tinf
use neutral_perturbations, only: clear_dneu
use gemini3d, only: fluidvar_pointers,electrovar_pointers,intvars

implicit none (type, external)
private
public :: init_procgrid, outdir_fullgridvaralloc, read_grid_in, get_initial_state, &
            BGfield_Lagrangian, check_dryrun, check_fileoutput,  &
            get_initial_drifts, init_Efieldinput, pot2perpfield_in, init_neutralperturb_in, dt_select, &
            neutral_atmos_wind_update, neutral_perturb, electrodynamics_in, check_finite_output_in, &
            halo_interface_vels_allspec_in, set_global_boundaries_allspec_in, halo_allparams_in, &
            RK2_prep_mpi_allspec_in,get_gavg_Tinf_in, clear_dneu_in, mpisetup, mpiparms

real(wp), dimension(:,:,:), allocatable :: Phiall    ! full grid potential, only root allocates
real(wp), parameter :: dtscale=2                     ! controls how rapidly the time step is allowed to change

contains
  !> call the mpi setup from our module
  subroutine mpisetup_in()
    call mpisetup()
  end subroutine mpisetup_in


  !> my id and the total number of processes for this run
  subroutine mpiparms(myid,lid)
    integer, intent(inout) :: myid,lid

    myid=mpi_cfg%myid
    lid=mpi_cfg%lid
  end subroutine mpiparms


  !> create output directory and allocate full grid potential storage
  subroutine outdir_fullgridvaralloc(cfg,lx1,lx2all,lx3all)
    type(gemini_cfg), intent(in) :: cfg
    integer, intent(in) :: lx1,lx2all,lx3all

    !> create a place, if necessary, for output datafiles
    if (mpi_cfg%myid==0) then
      call create_outdir(cfg)
      allocate(Phiall(-1:lx1+2,-1:lx2all+2,-1:lx3all+2))
    end if
  end subroutine outdir_fullgridvaralloc


  !> read in the grid and distribute to workers
  subroutine read_grid_in(cfg,x)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), pointer, intent(inout) :: x

    call read_grid(cfg%indatsize,cfg%indatgrid,cfg%flagperiodic, x)
    !! read in a previously generated grid from filenames listed in input file
  end subroutine read_grid_in


  !> load initial conditions
  subroutine get_initial_state(cfg,fluidvars,x,UTsec,ymd,tdur)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), intent(inout) :: fluidvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(inout) :: UTsec
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: tdur
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    integer, dimension(3) :: ymdtmp
    real(wp) :: UTsectmp,ttmp
    character(:), allocatable :: filetmp

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
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
      print*, 'Starting from beginning of simulation...'
      call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,cfg%indatfile,ns,vs1,Ts,Phi,Phiall)
    end if
  end subroutine get_initial_state


  !> check whether file output should be done and complete it
  subroutine check_fileoutput(cfg,fluidvars,electrovars,intvars,t,tout,tglowout,tmilestone,flagoutput,ymd,UTsec)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(in) :: intvars
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: tout,tglowout,tmilestone
    integer, intent(inout) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp) :: tstart,tfin
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
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
      call output_aur(cfg%outdir, cfg%flagglow, ymd, UTsec, intvars%iver, cfg%out_format)
      if (mpi_cfg%myid==0) then
        call cpu_time(tfin)
        print *, 'Auroral output done for time step:  ',t,' in cpu_time of: ',tfin-tstart
      end if
      tglowout = tglowout + cfg%dtglowout
    end if
  end subroutine


  !> check whether a dryrun simulation was done
  subroutine check_dryrun(cfg)
    type(gemini_cfg), intent(in) :: cfg
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


  !> prep simulation for use of Lagrangian grid, if needed
  subroutine BGfield_Lagrangian(x,v2grid,v3grid)
    class(curvmesh), intent(in) :: x
    real(wp), intent(inout) :: v2grid,v3grid
    real(wp), dimension(:,:,:), allocatable :: E01,E02,E03

    allocate(E01(lx1,lx2,lx3),E02(lx1,lx2,lx3),E03(lx1,lx2,lx3))    ! allocate without ghost cells
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
      E1(1:lx1,1:lx2,1:lx3) = E1(1:lx1,1:lx2,1:lx3) + E01    ! FIXME: this is before dist fields are computed???
      E2(1:lx1,1:lx2,1:lx3) = E2(1:lx1,1:lx2,1:lx3) + E02
      E3(1:lx1,1:lx2,1:lx3) = E3(1:lx1,1:lx2,1:lx3) + E03
    end if
    deallocate(E01,E02,E03)

    if (mpi_cfg%myid==0) then
      print*, 'Recomputed initial fields including background (BGfieldLagrangian):'
      print*, '    ',minval(E1(1:lx1,1:lx2,1:lx3)),maxval(E1(1:lx1,1:lx2,1:lx3))
      print*, '    ',minval(E2(1:lx1,1:lx2,1:lx3)),maxval(E2(1:lx1,1:lx2,1:lx3))
      print*, '    ',minval(E3(1:lx1,1:lx2,1:lx3)),maxval(E3(1:lx1,1:lx2,1:lx3))
    end if
  end subroutine BGfield_Lagrangian


  !> initial drifts at the start of the simulation
  subroutine get_initial_drifts(cfg,x,fluidvars,fluidauxvars,electrovars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    real(wp), dimension(:,:,:), allocatable :: sig0,sigP,sigH,sigPgrav,sigHgrav
    real(wp), dimension(:,:,:,:), allocatable :: muP,muH,nusn
    integer :: lx1,lx2,lx3,lsp
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    ! calculate drifts
    lx1=x%lx1; lx2=x%lx2; lx3=x%lx3; lsp=size(ns,4);
    allocate(sig0(lx1,lx2,lx3),sigP(lx1,lx2,lx3),sigH(lx1,lx2,lx3),sigPgrav(lx1,lx2,lx3),sigHgrav(lx1,lx2,lx3))
    allocate(muP(lx1,lx2,lx3,lsp),muH(lx1,lx2,lx3,lsp),nusn(lx1,lx2,lx3,lsp))
    call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
    call velocities(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)
    deallocate(sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
    if(mpi_cfg%myid==0) then
      print*, 'Recomputed initial drifts:  '
      print*, '    ',minval(vs2(1:lx1,1:lx2,1:lx3,1:lsp)),maxval(vs2(1:lx1,1:lx2,1:lx3,1:lsp))
      print*, '    ',minval(vs3(1:lx1,1:lx2,1:lx3,1:lsp)),maxval(vs3(1:lx1,1:lx2,1:lx3,1:lsp))
    end if
  end subroutine get_initial_drifts


  !> initialize the process gridf for this simulation
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


  !> initialize electric field input data
  subroutine init_Efieldinput_in(x,dt,t,intvars,ymd,UTsec)
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt,t
    type(gemini_work), intent(inout) :: intvars
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call init_Efieldinput(dt,t,cfg,ymd,UTsec,x,intvars%efield)
  end subroutine init_Efieldinput_in


  !> convert potential to electric field by differentiating
  subroutine pot2perpfield_in(x,electrovars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(inout) :: electrovars
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi

    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    E1 = 0
    call pot2perpfield(Phi,x,E2,E3)
    if(mpi_cfg%myid==0) then
      print '(A)', 'Recomputed initial dist. fields via pot2perpfield:'
      print*, '    gemini ',minval(E1(1:lx1,1:lx2,1:lx3)),maxval(E1(1:lx1,1:lx2,1:lx3))
      print*, '    gemini ',minval(E2(1:lx1,1:lx2,1:lx3)),maxval(E2(1:lx1,1:lx2,1:lx3))
      print*, '    gemini ',minval(E3(1:lx1,1:lx2,1:lx3)),maxval(E3(1:lx1,1:lx2,1:lx3))
      print*, '    gemini ',minval(Phi(1:lx1,1:lx2,1:lx3)),maxval(Phi(1:lx1,1:lx2,1:lx3))
    end if
  end subroutine pot2perpfield_in


  !> initialize neutral perturbations
  subroutine init_neutralperturb_in(dt,x,intvars,ymd,UTsec)
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    type(gemini_work), intent(inout) :: intvars
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call init_neutralperturb(cfg,x,dt,ymd,UTsec,intvars%atmosperturb)
  end subroutine init_neutralperturb_in


  !> select time step and throttle if changing too rapidly
  subroutine dt_select(cfg,x,fluidvars,fluidauxvars,it,t,tout,tglowout,dt)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: it
    real(wp), intent(in) :: t,tout,tglowout
    real(wp), intent(inout) :: dt
    real(wp) :: dtprev
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !> save prior time step
    dtprev = dt

    !> time step calculation, requires workers to report their most stringent local stability constraint
    call dt_comm(t,tout,tglowout,cfg,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x,dt)

    !> do not allow the time step to change too rapidly
    if (it>1) then
      if(dt/dtprev > dtscale) then
        !! throttle how quickly we allow dt to increase
        dt=dtscale*dtprev
        if (mpi_cfg%myid == 0) then
          print '(A,EN14.3)', 'Throttling dt to:  ',dt
        end if
      end if
    end if
  end subroutine dt_select


  !> apply neutral perturbations/background and assign to main code variables
  subroutine neutral_atmos_wind_update(intvars,v2grid,v3grid)
    type(gemini_work), intent(inout) :: intvars
    real(wp), intent(in) :: v2grid,v3grid

    call neutral_denstemp_update(intvars%atmos%nn,intvars%atmos%Tn)
    call neutral_wind_update(intvars%atmos%vn1,intvars%atmos%vn2,intvars%atmos%vn3,v2grid,v3grid)
  end subroutine neutral_atmos_wind_update


  !> compute neutral perturbations and apply to main code variables
  subroutine neutral_perturb_in(cfg,intvars,x,dt,t,ymd,UTsec,v2grid,v3grid)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: v2grid,v3grid

    call neutral_perturb(cfg,dt,cfg%dtneu,t,ymd,UTsec,x,v2grid,v3grid,intvars%atmos)
    call check_finite_pertub(cfg%outdir, t, mpi_cfg%myid, nn, Tn, vn1, vn2, vn3)
  end subroutine neutral_perturb_in


  !> call electrodynamics solution
  subroutine electrodynamics_in(fluidvars,fluidauxvars,electrovars,x,it,t,dt,ymd,UTsec)
    real(wp), dimension(:,:,:,:), intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), intent(in) :: electrovars
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: it
    real(wp), intent(in) :: t,dt
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    ! E&M solves
    call electrodynamics(it,t,dt,nn,vn2,vn3,Tn,cfg,ns,Ts,vs1,B1,vs2,vs3,x,E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)
  end subroutine electrodynamics_in


  !> check main state variables for finiteness
  subroutine check_finite_output_in(cfg,fluidvars,electrovars,t)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    real(wp), intent(in) :: t
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    call check_finite_output(cfg%outdir,t,mpi_cfg%myid,vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)
  end subroutine check_finite_output_in


  !> haloing for computing cell interface velocities
  subroutine halo_interface_vels_allspec_in(x,fluidvars,lsp)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(inout) :: fluidvars    
    integer, intent(in) :: lsp
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call halo_interface_vels_allspec(x%flagper,vs2,vs3,lsp)
  end subroutine halo_interface_vels_allspec_in


  !> enforce global boundary conditions
  subroutine set_global_boundaries_allspec_in(x,fluidvars,fluidauxvars,intvars,lsp)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(inout) :: fluidvars    
    real(wp), dimension(:,:,:,:), intent(inout) :: fluidauxvars
    integer, intent(in) :: lsp
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    ! fix global boundaries, as needed
    call set_global_boundaries_allspec(x%flagper,ns,rhovs1,vs1,vs2,vs3,rhoes,vs1i,lsp)
  end subroutine set_global_boundaries_allspec_in


  !> halo all advected parameters
  subroutine halo_allparams_in(x,fluidvars,fluidauxvars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(inout) :: fluidvars    
    real(wp), dimension(:,:,:,:), intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    ! halo all fluid parameters
    call halo_allparams(ns,rhovs1,rhoes,x%flagper)
  end subroutine halo_allparams_in


  !> prepare/halo data for compression substep
  subroutine RK2_prep_mpi_allspec_in(x,fluidvars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(inout) :: fluidvars    
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call RK2_prep_mpi_allspec(vs1,vs2,vs3,x%flagper)
  end subroutine RK2_prep_mpi_allspec_in


  !> agree on average value of gravity and exospheric temp
  subroutine get_gavg_Tinf_in(gavg,Tninf)
    real(wp), intent(inout) :: gavg,Tninf

    call get_gavg_Tinf(gavg,Tninf)
  end subroutine get_gavg_Tinf_in


  !> deallocate module storage for neutral perturbations
  subroutine clear_dneu_in(intvars)
    call clear_dneu()
  end subroutine clear_dneu_in
end module gemini3d_mpi

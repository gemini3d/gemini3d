! Copyright 2021 Matthew Zettergren

! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

program Gemini3D_main
!! a main program illustrating use of gemini library to conduct an ionospheric simulation
use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_bool, c_float
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use sanity_check, only : check_finite_output, check_finite_pertub
use phys_consts, only : lnchem, lwave, lsp, wp, debug
use grid, only: lx1,lx2,lx3,lx2all,lx3all
use grid_mpi, only: read_grid,grid_drift
use meshobj, only: curvmesh
use config, only : gemini_cfg
use io, only : input_plasma,create_outdir,create_outdir_aur
use mpimod, only : mpisetup, mpibreakdown, mpi_manualgrid, process_grid_auto, mpi_cfg
use multifluid, only : sweep3_allparams,sweep1_allparams,sweep2_allparams,source_loss_allparams,VNRicht_artvisc,compression, &
            energy_diffusion,impact_ionization,clean_param,rhoe2T,T2rhoe,rhov12v1,v12rhov1
use ionization_mpi, only: get_gavg_Tinf
use multifluid_mpi, only: halo_allparams
use msis_interface, only : msisinit
use neutral, only : neutral_atmos,make_neuBG,init_neutralBG,neutral_winds,clear_neuBG
use neutral_perturbations, only: init_neutralperturb,neutral_perturb,clear_dneu,neutral_denstemp_update,neutral_wind_update
use potentialBCs_mumps, only: init_Efieldinput
use potential_comm,only : electrodynamics,pot2perpfield
use precipBCs_mod, only: init_precipinput
use temporal, only : dt_comm
use timeutils, only: dateinc, find_lastdate
use advec, only: interface_vels_allspec
use advec_mpi, only: halo_interface_vels_allspec,set_global_boundaries_allspec
use sources_mpi, only: RK2_prep_mpi_allspec
use gemini3d, only: c_params,cli_config_gridsize,gemini_alloc,gemini_dealloc,cfg,x
use gemini3d_mpi, only: init_procgrid,outdir_fullgridvaralloc,get_initial_state,BGfield_Lagrangian, &
                          check_dryrun,check_fileoutput,get_initial_drifts

implicit none (type, external)
external :: mpi_init

integer(c_int) :: lid2in, lid3in
character(8) :: date
character(10) :: time
integer :: ierr
type(c_params) :: p

!> initialize mpi
call mpi_init(ierr)
if (ierr/=0) error stop 'gemini.bin: failed mpi_init'
p%fortran_cli = .true.
p%out_dir(1) = c_null_char
lid2in = -1
lid3in = -1

!! out_dir, lid2in, lid3in, are ignored when fortran_cli=.true.
call gemini_main(p, lid2in, lid3in)

!> shut down mpi
ierr = mpibreakdown()

if (ierr /= 0) then
  write(stderr, *) 'GEMINI: abnormal MPI shutdown code', ierr, 'Process #', mpi_cfg%myid,' /',mpi_cfg%lid-1
  error stop
endif

call date_and_time(date,time)
print '(/,A,I0,A,I0,A)', 'GEMINI normal termination, Process # ', mpi_cfg%myid,' / ',mpi_cfg%lid-1, ' at ' // date // 'T' // time

contains
  subroutine gemini_main(p, lid2in, lid3in)  bind(C)
    !! NOTE: if use_cli=.true., then {out_dir, lid2in, lid3in} are ignored and CLI is used instead.
    type(c_params), intent(in) :: p
    !! output directory for Gemini3D to write simulation data to (can be large files GB, TB, ...)
    integer(c_int), intent(inout) :: lid2in, lid3in  !< inout to allow optional CLI
    integer :: ierr
    logical :: exists
    !> VARIABLES READ IN FROM CONFIG FILE
    real(wp) :: UTsec
    !! UT (s)
    integer, dimension(3) :: ymd
    !! year, month, day (current, not to be confused with starting year month and day in gemini_cfg structure)
      !> STATE VARIABLES
    !> MZ note:  it is likely that there could be a plasma and neutral derived type containing these data...  May be worth considering in a refactor...
    real(wp), dimension(:,:,:,:), allocatable :: ns,vs1,vs2,vs3,Ts
    !! fluid state variables
    real(wp), dimension(:,:,:), allocatable :: E1,E2,E3,J1,J2,J3,Phi
    !! electrodynamic state variables
    real(wp), dimension(:,:,:), allocatable :: rhov2,rhov3,B1,B2,B3
    !! inductive state vars. (for future use - except for B1 which is used for the background field)
    real(wp), dimension(:,:,:), allocatable :: rhom,v1,v2,v3
    !! inductive auxiliary
    real(wp), dimension(:,:,:,:), allocatable :: nn
    !! neutral density array
    real(wp), dimension(:,:,:), allocatable :: Tn,vn1,vn2,vn3
    !! neutral temperature and velocities
    real(wp), dimension(:,:,:), allocatable :: Phiall
    !! full-grid potential solution.  To store previous time step value
    real(wp), dimension(:,:,:), allocatable :: iver
    !! integrated volume emission rate of aurora calculated by GLOW
    !TEMPORAL VARIABLES
    real(wp) :: t=0, dt=1e-6_wp,dtprev
    !! time from beginning of simulation (s) and time step (s)
    real(wp) :: tout
    !! time for next output and time between outputs
    real(wp) :: tstart,tfin
    !! temp. vars. for measuring performance of code blocks
    integer :: it,isp,iupdate
    !! time and species loop indices
    real(wp) :: tneuBG !for testing whether we should re-evaluate neutral background
      !> WORK ARRAYS
    real(wp), allocatable :: dl1,dl2,dl3     !these are grid distances in [m] used to compute Courant numbers
    real(wp) :: tglowout,tdur
    !! time for next GLOW output
    !> TO CONTROL THROTTLING OF TIME STEP
    real(wp), parameter :: dtscale=2
    !> Temporary variable for toggling full vs. other output
    integer :: flagoutput
    real(wp) :: tmilestone = 0
    !> Describing Lagrangian grid (if used)
    real(wp) :: v2grid,v3grid
    character(*), parameter :: msis2_param_file = "msis20.parm"
    
    !> initialize message passing
    call mpisetup(); if(mpi_cfg%lid < 1) error stop 'number of MPI processes must be >= 1. Was MPI initialized properly?'
    call cli_config_gridsize(p,lid2in,lid3in)
    
    !> MPI gridding cannot be done until we know the grid size, and needs to be done before we distribute pieces of the grid
    !    to workers
    call init_procgrid(lx2all,lx3all,lid2in,lid3in)
    
    !> load the grid data from the input file
    call read_grid(cfg%indatsize,cfg%indatgrid,cfg%flagperiodic, x)
    !! read in a previously generated grid from filenames listed in input file
    
    !> Allocate space for solutions
    call gemini_alloc(ns,vs1,vs2,vs3,Ts,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom, &
                        E1,E2,E3,J1,J2,J3,Phi,nn,Tn,vn1,vn2,vn3,iver)
 
    !> root creates a place to put output and allocates any needed fullgrid arrays for plasma state variables
    call outdir_fullgridvaralloc(Phiall,lx1,lx2all,lx3all)
    
    !> Set initial time variables to simulation; this requires detecting whether we are trying to restart a simulation run
    call get_initial_state(ns,vs1,Ts,Phi,Phiall,UTsec,ymd,tdur)

    !> Initialize some variables need for time stepping and output 
    it = 1; t = 0; tout = t; tglowout = t; tneuBG=t
    
    !ROOT/WORKERS WILL ASSUME THAT THE MAGNETIC FIELDS AND PERP FLOWS START AT ZERO
    !THIS KEEPS US FROM HAVING TO HAVE FULL-GRID ARRAYS FOR THESE STATE VARS (EXCEPT
    !FOR IN OUTPUT FNS.).  IF A SIMULATIONS IS DONE WITH INERTIAL CAPACITANCE THERE
    !WILL BE A FINITE AMOUNT OF TIME FOR THE FLOWS TO 'START UP', BUT THIS SHOULDN'T
    !BE TOO MUCH OF AN ISSUE.  WE ALSO NEED TO SET THE BACKGROUND MAGNETIC FIELD STATE
    !VARIABLE HERE TO WHATEVER IS SPECIFIED IN THE GRID STRUCTURE (THESE MUST BE CONSISTENT)
    rhov2 = 0; rhov3 = 0; v2 = 0; v3 = 0; B2 = 0; B3 = 0; B1(1:lx1,1:lx2,1:lx3) = x%Bmag
    !! this assumes that the grid is defined s.t. the x1 direction corresponds
    !! to the magnetic field direction (hence zero B2 and B3).
    
    !> Inialize neutral atmosphere, note the use of fortran's weird scoping rules to avoid input args.  Must occur after initial time info setup
    if(mpi_cfg%myid==0) print*, 'Priming electric field input'
    call init_Efieldinput(dt,t,cfg,ymd,UTsec,x)
    
    !> Recompute electrodynamic quantities needed for restarting
    !> these do not include background
    E1 = 0
    call pot2perpfield(Phi,x,E2,E3)    
    if(mpi_cfg%myid==0) then
      print '(A)', 'Recomputed initial dist. fields:'
      print*, '    gemini ',minval(E1),maxval(E1)
      print*, '    gemini ',minval(E2),maxval(E2)
      print*, '    gemini ',minval(E3),maxval(E3)
    end if
    !> Get the background electric fields and compute the grid drift speed if user selected lagrangian grid, add to total field
    call BGfield_Lagrangian(v2grid,v3grid,E1,E2,E3) 
    if (mpi_cfg%myid==0) then    
      print*, 'Recomputed initial fields including background:'
      print*, '    ',minval(E1),maxval(E1)
      print*, '    ',minval(E2),maxval(E2)
      print*, '    ',minval(E3),maxval(E3)
    end if
    
    !> Precipitation input setup
    if(mpi_cfg%myid==0) print*, 'Priming precipitation input'
    call init_precipinput(dt,t,cfg,ymd,UTsec,x)
    
    !> Neutral atmosphere setup
    if(cfg%msis_version == 20) then
      inquire(file=msis2_param_file, exist=exists)
      if(.not.exists) error stop 'could not find MSIS 2 parameter file ' // msis2_param_file // &
        ' this file must be in the same directory as gemini.bin, and run from that directory. ' // &
        'This limitation comes from how MSIS 2.x is coded internally.'
      call msisinit(parmfile=msis2_param_file)
    end if
    if(mpi_cfg%myid==0) print*, 'Computing background and priming neutral perturbation input (if used)'
    call init_neutralBG(dt,t,cfg,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
    call init_neutralperturb(cfg,x,dt,ymd,UTsec)

 
    !> Recompute drifts and make some decisions about whether to invoke a Lagrangian grid
    call get_initial_drifts(nn,Tn,vn1,vn2,vn3,ns,Ts,vs1,vs2,vs3,B1,E2,E3)
    if(mpi_cfg%myid==0) then
      print*, 'Recomputed initial drifts:  '
      print*, '    ',minval(vs2(1:lx1,1:lx2,1:lx3,1:lsp)),maxval(vs2(1:lx1,1:lx2,1:lx3,1:lsp))
      print*, '    ',minval(vs3(1:lx1,1:lx2,1:lx3,1:lsp)),maxval(vs3(1:lx1,1:lx2,1:lx3,1:lsp))
    end if
    
    !> control update rate from excessive console printing
    !! considering small vs. large simulations
    !! these are arbitrary levels, so feel free to finesse
    if (lx1*lx2*lx3 < 20000) then
      iupdate = 50
    elseif (lx1*lx2*lx3 < 100000) then
      iupdate = 10
    else
      iupdate = 1
    endif
    
    !> Main time loop
    main : do while (t < tdur)
      !> time step calculation, requires workers to report their most stringent local stability constraint
      dtprev = dt
      call dt_comm(t,tout,tglowout,cfg,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x,dt)
      if (it>1) then
        if(dt/dtprev > dtscale) then
          !! throttle how quickly we allow dt to increase
          dt=dtscale*dtprev
          if (mpi_cfg%myid == 0) then
            print '(A,EN14.3)', 'Throttling dt to:  ',dt
          end if
        end if
      end if
    
      !> get neutral background
      if ( it/=1 .and. cfg%flagneuBG .and. t>tneuBG) then     !we dont' throttle for tneuBG so we have to do things this way to not skip over...
        call cpu_time(tstart)
        call neutral_atmos(ymd,UTsec,x%glat,x%glon,x%alt,cfg%activ,cfg%msis_version)
        call neutral_winds(ymd, UTsec, Ap=cfg%activ(3), x=x)
        call neutral_denstemp_update(nn,Tn)    ! empirical model calls don't automatically assign state to main variables
        call neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid)
        tneuBG=tneuBG+cfg%dtneuBG
        if (mpi_cfg%myid==0) then
          call cpu_time(tfin)
          print *, 'Neutral background at time:  ',t,' calculated in time:  ',tfin-tstart
        end if
      end if
    
      !> get neutral perturbations
      if (cfg%flagdneu==1) then
        call cpu_time(tstart)
        call neutral_perturb(cfg,dt,cfg%dtneu,t,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
        call check_finite_pertub(cfg%outdir, t, mpi_cfg%myid, nn, Tn, vn1, vn2, vn3)
        if (mpi_cfg%myid==0 .and. debug) then
          call cpu_time(tfin)
          print *, 'Neutral perturbations calculated in time:  ',tfin-tstart
        endif
      end if
    
      !> compute potential solution
      call cpu_time(tstart)
      call electrodynamics(it,t,dt,nn,vn2,vn3,Tn,cfg,ns,Ts,vs1,B1,vs2,vs3,x,E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)
      if (mpi_cfg%myid==0 .and. debug) then
        call cpu_time(tfin)
        print *, 'Electrodynamics total solve time:  ',tfin-tstart
      endif
    
      !> update fluid variables
      if (mpi_cfg%myid==0 .and. debug) call cpu_time(tstart)
      call fluid_adv(ns,vs1,Ts,vs2,vs3,J1,E1,cfg,t,dt,x,nn,vn1,vn2,vn3,Tn,iver,ymd,UTsec, first=(it==1) )
      if (mpi_cfg%myid==0 .and. debug) then
        call cpu_time(tfin)
        print *, 'Multifluid total solve time:  ',tfin-tstart
      endif
    
      !> Sanity check key variables before advancing
      ! FIXME: for whatever reason, it is just a fact that vs1 has trash in ghost cells after fluid_adv; I don't know why...
      call check_finite_output(cfg%outdir, t, mpi_cfg%myid, vs2,vs3,ns,vs1,Ts, Phi,J1,J2,J3)
    
      !> update time variables
      it = it + 1
      t = t + dt
      if (mpi_cfg%myid==0 .and. debug) print *, 'Moving on to time step (in sec):  ',t,'; end time of simulation:  ',cfg%tdur
      call dateinc(dt,ymd,UTsec)
      if (mpi_cfg%myid==0 .and. (modulo(it, iupdate) == 0 .or. debug)) then
        !! print every 10th time step to avoid extreme amounts of console printing
        print '(A,I4,A1,I0.2,A1,I0.2,A1,F12.6,A5,F8.6)', 'Current time ',ymd(1),'-',ymd(2),'-',ymd(3),' ',UTsec,'; dt=',dt
      endif
    
      !> see if we are doing a dry run and exit program if so
      call check_dryrun()
    
      !> File output
      call check_fileoutput(t,tout,tglowout,tmilestone,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3,iver)
    end do main
    
    !> deallocate variables and module data
    call gemini_dealloc(ns,vs1,vs2,vs3,Ts,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom,E1,E2,E3,J1,J2,J3,Phi,nn,Tn,vn1,vn2,vn3,iver)
    if (mpi_cfg%myid==0) deallocate(Phiall)
    call clear_neuBG()
    call clear_dneu()
  end subroutine gemini_main
  
  
  !> this advances the fluid soluation by time interval dt
  subroutine fluid_adv(ns,vs1,Ts,vs2,vs3,J1,E1,cfg,t,dt,x,nn,vn1,vn2,vn3,Tn,iver,ymd,UTsec,first)
    !! J1 needed for heat conduction; E1 for momentum equation
    !! THIS SUBROUTINE ADVANCES ALL OF THE FLUID VARIABLES BY TIME STEP DT.
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  ns,vs1,Ts
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  vs2,vs3
    real(wp), dimension(:,:,:), intent(in) :: J1
    !! needed for thermal conduction in electron population
    real(wp), dimension(:,:,:), intent(inout) :: E1
    !! will have ambipolar field added into it in this procedure...
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: t,dt
    class(curvmesh), intent(in) :: x
    !! grid structure variable
    real(wp), dimension(:,:,:,:), intent(in) :: nn
    real(wp), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    logical, intent(in) :: first  !< first time step
    real(wp), dimension(:,:,:), intent(inout) :: iver
    !! intent(out)
    integer :: isp
    real(wp) :: tstart,tfin
    real(wp) :: f107,f107a
    real(wp), dimension(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2,size(ns,4)) ::  rhovs1,rhoes
    real(wp), dimension(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2) :: param
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: paramtrim
    real(wp), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4,size(ns,4)) :: vs1i
    real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-3,1:size(vs1,3)-4,size(ns,4)) :: vs2i
    real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-3,size(ns,4)) :: vs3i
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)) :: Q    ! artificial viscosity
    real(wp) :: gavg,Tninf
    
    ! cfg arrays can be confusing, particularly f107, so assign to sensible variable name here
    f107=cfg%activ(2)
    f107a=cfg%activ(1)
    
    ! Prior to advection substep convert velocity and temperature to momentum and enegy density (which are local to this procedure)
    call v12rhov1(ns,vs1,rhovs1)
    call T2rhoe(ns,Ts,rhoes) 
   
    ! advection substep for all species
    call cpu_time(tstart)
    call halo_interface_vels_allspec(x%flagper,vs2,vs3,vs2i,vs3i,lsp)
    call interface_vels_allspec(vs1,vs2,vs3,vs1i,vs2i,vs3i,lsp)    ! needs to happen regardless of ions v. electron due to energy eqn.
    call set_global_boundaries_allspec(x%flagper,ns,rhovs1,vs1,vs2,vs3,rhoes,vs1i,lsp)
    call halo_allparams(ns,rhovs1,rhoes,x%flagper)
    call sweep3_allparams(dt,x,vs3i,ns,rhovs1,rhoes)
    call sweep1_allparams(dt,x,vs1i,ns,rhovs1,rhoes)
    call halo_allparams(ns,rhovs1,rhoes,x%flagper)
    call sweep2_allparams(dt,x,vs2i,ns,rhovs1,rhoes)
    call rhov12v1(ns,rhovs1,vs1)
    call cpu_time(tfin)
    if (mpi_cfg%myid==0 .and. debug) then
      print *, 'Completed advection substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if
    
    ! post advection filling of null cells
    call clean_param(x,1,ns)
    call clean_param(x,2,vs1)
    
    ! Compute artifical viscosity and then execute compression calculation
    call cpu_time(tstart)
    call VNRicht_artvisc(ns,vs1,Q)
    call RK2_prep_mpi_allspec(vs1,vs2,vs3,x%flagper)
    call compression(dt,x,vs1,vs2,vs3,Q,rhoes)   ! this applies compression substep and then converts back to temperature
    call rhoe2T(ns,rhoes,Ts)
    call clean_param(x,3,Ts)
    call cpu_time(tfin)
    if (mpi_cfg%myid==0 .and. debug) then
      print *, 'Completed compression substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if
    
    ! Energy diffusion (thermal conduction) substep
    call cpu_time(tstart)
    call energy_diffusion(dt,x,ns,Ts,J1,nn,Tn,cfg%diffsolvetype,cfg%Teinf)
    call cpu_time(tfin)
    if (mpi_cfg%myid==0 .and. debug) then
      print *, 'Completed energy diffusion substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if
    
    ! cleanup and convert to specific internal energy density for sources substeps
    call clean_param(x,3,Ts)
    call T2rhoe(ns,Ts,rhoes)

    !> all workers need to "agree" on a gravity and exospheric temperature
    call get_gavg_Tinf(gavg,Tninf)

    !> solve all source/loss processes
    call source_loss_allparams(dt,t,cfg,ymd,UTsec,x,E1,Q,f107a,f107,nn,vn1,vn2,vn3, &
                                     Tn,first,ns,rhovs1,rhoes,vs1,vs2,vs3,Ts,iver,gavg,Tninf)
  
    ! density to be cleaned after source/loss
    call clean_param(x,3,Ts)
    call clean_param(x,2,vs1)
    call clean_param(x,1,ns)
    
    !should the electron velocity be recomputed here now that densities have changed...
  end subroutine fluid_adv
end program

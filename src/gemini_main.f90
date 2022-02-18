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
use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_bool, c_float, c_ptr
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use phys_consts, only : lnchem, lwave, lsp, wp, debug
!use grid, only: lx1,lx2,lx3,lx2all,lx3all
use mpimod, only : mpisetup, mpibreakdown, mpi_manualgrid, process_grid_auto, mpi_cfg

!> main gemini libraries
use gemini3d, only: c_params,cli_config_gridsize,gemini_alloc,gemini_dealloc,cfg,x,init_precipinput_C,msisinit_C, &
                      set_start_values, init_neutralBG_C, set_update_cadence, neutral_atmos_winds_C, get_solar_indices_C, &
                      v12rhov1_C,T2rhoe_C,interface_vels_allspec_C, sweep3_allparams_C, sweep1_allparams_C, sweep2_allparams_C, &
                      rhov12v1_C, VNRicht_artvisc_C, compression_C, rhoe2T_C, clean_param_C, energy_diffusion_C, &
                      source_loss_allparams_C,clear_neuBG_C,dateinc_C,get_subgrid_size_C, get_fullgrid_size_C
use gemini3d_mpi, only: init_procgrid,outdir_fullgridvaralloc,read_grid_C,get_initial_state,BGfield_Lagrangian, &
                          check_dryrun,check_fileoutput,get_initial_drifts,init_Efieldinput_C,pot2perpfield_C, &
                          init_neutralperturb_C, dt_select_C, neutral_atmos_wind_update_C, neutral_perturb_C, &
                          electrodynamics_C, check_finite_output_C, halo_interface_vels_allspec_C, &
                          set_global_boundaries_allspec_C, halo_allparams_C, RK2_prep_mpi_allspec_C, get_gavg_Tinf_C, &
                          clear_dneu_C

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
    !> VARIABLES READ IN FROM CONFIG FILE
    real(wp) :: UTsec
    !! UT (s)
    integer, dimension(3) :: ymd
    !! year, month, day (current, not to be confused with starting year month and day in gemini_cfg structure)
      !> STATE VARIABLES
    !> MZ note:  it is likely that there could be a plasma and neutral derived type containing these data...  May be worth considering in a refactor...
    type(c_ptr) :: fluidvarsC    ! large data block for holding gemini fluid variables
    type(c_ptr) :: fluidauxvarsC   ! aux variables memory block for fluid parms
    type(c_ptr) :: electrovarsC   ! large block of data for electrodynamic variables (excluding mag. fields)
    !TEMPORAL VARIABLES
    real(wp) :: t=0, dt=1e-6_wp
    !! time from beginning of simulation (s) and time step (s)
    real(wp) :: tout
    !! time for next output and time between outputs
    real(wp) :: tstart,tfin
    !! temp. vars. for measuring performance of code blocks
    integer :: it,iupdate
    !! time and species loop indices
    real(wp) :: tneuBG !for testing whether we should re-evaluate neutral background
      !> WORK ARRAYS
    real(wp) :: tglowout,tdur
    !! time for next GLOW output
    !> Temporary variable for toggling full vs. other output
    integer :: flagoutput
    real(wp) :: tmilestone = 0
    !> Describing Lagrangian grid (if used)
    real(wp) :: v2grid,v3grid
    integer :: lx1,lx2,lx3,lx2all,lx3all
    
    !> initialize message passing
    call mpisetup(); if(mpi_cfg%lid < 1) error stop 'number of MPI processes must be >= 1. Was MPI initialized properly?'
    call cli_config_gridsize(p,lid2in,lid3in)
    call get_fullgrid_size_C(lx1,lx2all,lx3all)
    
    !> MPI gridding cannot be done until we know the grid size, and needs to be done before we distribute pieces of the grid
    !    to workers
    call init_procgrid(lx2all,lx3all,lid2in,lid3in)
    
    !> load the grid data from the input file and store in gemini module
    call read_grid_C()
    call get_subgrid_size_C(lx1,lx2,lx3)
    
    !> Allocate space for solutions
    call gemini_alloc(fluidvarsC,fluidauxvarsC,electrovarsC)

    !> root creates a place to put output and allocates any needed fullgrid arrays for plasma state variables
    call outdir_fullgridvaralloc(lx1,lx2all,lx3all)
    
    !> Set initial time variables to simulation; this requires detecting whether we are trying to restart a simulation run
    call get_initial_state(UTsec,ymd,tdur)

    !> initialize time stepping and some aux variables
    call set_start_values(it,t,tout,tglowout,tneuBG)
    
    !> Inialize neutral atmosphere, note the use of fortran's weird scoping rules to avoid input args.  Must occur after initial time info setup
    if(mpi_cfg%myid==0) print*, 'Priming electric field input'
    call init_Efieldinput_C(dt,t,ymd,UTsec)
    
    !> Recompute electrodynamic quantities needed for restarting
    !> these do not include background
    call pot2perpfield_C()

    !> Get the background electric fields and compute the grid drift speed if user selected lagrangian grid, add to total field
    call BGfield_Lagrangian(v2grid,v3grid) 
    
    !> Precipitation input setup
    if(mpi_cfg%myid==0) print*, 'Priming precipitation input'
    call init_precipinput_C(dt,t,ymd,UTsec)
    
    !> Neutral atmosphere setup
    if(mpi_cfg%myid==0) print*, 'Computing background and priming neutral perturbation input (if used)'
    call msisinit_C()
    call init_neutralBG_C(dt,t,ymd,UTsec,v2grid,v3grid)
    call init_neutralperturb_C(dt,ymd,UTsec)

    !> Recompute drifts and make some decisions about whether to invoke a Lagrangian grid
    call get_initial_drifts()
    
    !> control rate of console printing
    call set_update_cadence(iupdate)
    
    !> Main time loop
    main : do while (t < tdur)
      call dt_select_C(it,t,tout,tglowout,dt)
    
      !> get neutral background
      if ( it/=1 .and. cfg%flagneuBG .and. t>tneuBG) then     !we dont' throttle for tneuBG so we have to do things this way to not skip over...
        call cpu_time(tstart)
        call neutral_atmos_winds_C(ymd,UTsec)   ! load background states into module variables
        call neutral_atmos_wind_update_C(v2grid,v3grid)    ! apply to variables in this program unit
        tneuBG=tneuBG+cfg%dtneuBG
        if (mpi_cfg%myid==0) then
          call cpu_time(tfin)
          print *, 'Neutral background at time:  ',t,' calculated in time:  ',tfin-tstart
        end if
      end if
    
      !> get neutral perturbations
      if (cfg%flagdneu==1) then
        call cpu_time(tstart)
        call neutral_perturb_C(dt,t,ymd,UTsec,v2grid,v3grid)
        if (mpi_cfg%myid==0 .and. debug) then
          call cpu_time(tfin)
          print *, 'Neutral perturbations calculated in time:  ',tfin-tstart
        endif
      end if
    
      !> compute potential solution
      call cpu_time(tstart)
      call electrodynamics_C(it,t,dt,ymd,UTsec)
      if (mpi_cfg%myid==0 .and. debug) then
        call cpu_time(tfin)
        print *, 'Electrodynamics total solve time:  ',tfin-tstart
      endif
    
      !> update fluid variables
      if (mpi_cfg%myid==0 .and. debug) call cpu_time(tstart)
      call fluid_adv(t,dt,ymd,UTsec, first=(it==1) )
      if (mpi_cfg%myid==0 .and. debug) then
        call cpu_time(tfin)
        print *, 'Multifluid total solve time:  ',tfin-tstart
      endif
    
      !> Sanity check key variables before advancing
      ! FIXME: for whatever reason, it is just a fact that vs1 has trash in ghost cells after fluid_adv; I don't know why...
      call check_finite_output_C(t)
    
      !> update time variables
      it = it + 1
      t = t + dt
      if (mpi_cfg%myid==0 .and. debug) print *, 'Moving on to time step (in sec):  ',t,'; end time of simulation:  ',cfg%tdur
      call dateinc_C(dt,ymd,UTsec)
      if (mpi_cfg%myid==0 .and. (modulo(it, iupdate) == 0 .or. debug)) then
        !! print every 10th time step to avoid extreme amounts of console printing
        print '(A,I4,A1,I0.2,A1,I0.2,A1,F12.6,A5,F8.6)', 'Current time ',ymd(1),'-',ymd(2),'-',ymd(3),' ',UTsec,'; dt=',dt
      endif
    
      !> see if we are doing a dry run and exit program if so
      call check_dryrun()
    
      !> File output
      call check_fileoutput(t,tout,tglowout,tmilestone,flagoutput,ymd,UTsec)
    end do main
    
    !> deallocate variables and module data
    call gemini_dealloc(fluidvarsC,fluidauxvarsC,electrovarsC)
    call clear_neuBG_C()
    call clear_dneu_C()
  end subroutine gemini_main
  
  
  !> this advances the fluid soluation by time interval dt
  subroutine fluid_adv(t,dt,ymd,UTsec,first)
    !! J1 needed for heat conduction; E1 for momentum equation
    !! THIS SUBROUTINE ADVANCES ALL OF THE FLUID VARIABLES BY TIME STEP DT.
    real(wp), intent(in) :: t,dt
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    logical, intent(in) :: first  !< first time step
    real(wp) :: tstart,tfin
    real(wp) :: f107,f107a
    real(wp) :: gavg,Tninf
    
    ! pull solar indices from module type
    call get_solar_indices_C(f107,f107a)

    ! Prior to advection substep convert velocity and temperature to momentum and enegy density (which are local to this procedure)
    call v12rhov1_C()
    call T2rhoe_C() 
   
    ! advection substep for all species
    call cpu_time(tstart)
    call halo_interface_vels_allspec_C(lsp)
    call interface_vels_allspec_C(lsp)    ! needs to happen regardless of ions v. electron due to energy eqn.
    call set_global_boundaries_allspec_C(lsp)
    call halo_allparams_C()
    call sweep3_allparams_C(dt)
    call sweep1_allparams_C(dt)
    call halo_allparams_C()
    call sweep2_allparams_C(dt)
    call rhov12v1_C()
    call cpu_time(tfin)
    if (mpi_cfg%myid==0 .and. debug) then
      print *, 'Completed advection substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if
    
    ! post advection filling of null cells
    call clean_param_C(1)
    call clean_param_C(2)
    
    ! Compute artifical viscosity and then execute compression calculation
    call cpu_time(tstart)
    call VNRicht_artvisc_C()
    call RK2_prep_mpi_allspec_C()     ! halos velocity so we can take a divergence without artifacts
    call compression_C(dt)   ! this applies compression substep and then converts back to temperature
    call rhoe2T_C()
    call clean_param_C(3)
    call cpu_time(tfin)
    if (mpi_cfg%myid==0 .and. debug) then
      print *, 'Completed compression substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if
    
    ! Energy diffusion (thermal conduction) substep
    call cpu_time(tstart)
    call energy_diffusion_C(dt)
    call cpu_time(tfin)
    if (mpi_cfg%myid==0 .and. debug) then
      print *, 'Completed energy diffusion substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if
    
    ! cleanup and convert to specific internal energy density for sources substeps
    call clean_param_C(3)
    call T2rhoe_C()

    !> all workers need to "agree" on a gravity and exospheric temperature
    call get_gavg_Tinf_C(gavg,Tninf)

    !> solve all source/loss processes
    call source_loss_allparams_C(dt,t,ymd,UTsec,f107a,f107,first,gavg,Tninf)
  
    ! density to be cleaned after source/loss
    call clean_param_C(3)
    call clean_param_C(2)
    call clean_param_C(1)
    
    !should the electron velocity be recomputed here now that densities have changed...
  end subroutine fluid_adv
end program

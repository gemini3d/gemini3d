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
use phys_consts, only : wp, debug
use mpi_f08, only: MPI_COMM_WORLD, mpi_init,mpi_finalize,mpi_comm_rank

!> type definitions
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg

!> main gemini libraries
use gemini3d, only: c_params,gemini_alloc,gemini_dealloc,init_precipinput_in,msisinit_in, &
                      set_start_values_auxtimevars, set_start_values_auxvars, init_neutralBG_in, &
                      set_update_cadence, neutral_atmos_winds, get_solar_indices, &
                      v12rhov1_in,T2rhoe_in,interface_vels_allspec_in, sweep3_allparams_in, &
                      sweep1_allparams_in, sweep2_allparams_in, &
                      rhov12v1_in, VNRicht_artvisc_in, compression_in, rhoe2T_in, clean_param_in, energy_diffusion_in, &
                      source_loss_allparams_in,dateinc_in,get_subgrid_size, get_fullgrid_size, &
                      get_config_vars, get_species_size, gemini_work, gemini_cfg_alloc, cli_in, read_config_in, &
                      gemini_cfg_dealloc, grid_size_in, gemini_double_alloc, gemini_work_alloc, gemini_double_dealloc, &
                      gemini_work_dealloc, set_global_boundaries_allspec_in
use gemini3d_mpi, only: init_procgrid,outdir_fullgridvaralloc,read_grid_in,get_initial_state,BGfield_Lagrangian, &
                          check_dryrun,check_fileoutput,get_initial_drifts,init_Efieldinput_in,pot2perpfield_in, &
                          init_neutralperturb_in, dt_select, neutral_atmos_wind_update, neutral_perturb_in, &
                          electrodynamics_in, check_finite_output_in, halo_interface_vels_allspec_in, &
                          halo_allparams_in, RK2_prep_mpi_allspec_in, get_gavg_Tinf_in, &
                          clear_dneu_in,mpisetup_in,mpiparms, calc_subgrid_size_in, halo_fluidvars_in, &
                          RK2_global_boundary_allspec_in

implicit none (type, external)

integer(c_int) :: lid2in, lid3in
character(8) :: date
character(10) :: time
integer :: ierr
type(c_params) :: p
integer :: myid

!> initialize mpi
call mpi_init()
p%fortran_cli = .true.
p%fortran_nml = .true.
p%out_dir(1) = c_null_char
lid2in = -1
lid3in = -1

!! out_dir, lid2in, lid3in, are ignored when fortran_cli=.true.
call gemini_main(p, lid2in, lid3in)

!> shut down mpi
call mpi_finalize(ierr)

if (ierr /= 0) then
  write(stderr, *) 'GEMINI: abnormal MPI shutdown code', ierr, 'Process #', myid
  error stop
endif

call date_and_time(date,time)
print '(/,A,I0,A,I0,A)', 'GEMINI normal termination, Process # ', myid,' at ' // date // 'T' // time

contains
  subroutine gemini_main(p, lid2in, lid3in)  bind(C)
    !! NOTE: if use_cli=.true., then {out_dir, lid2in, lid3in} are ignored and CLI is used instead.
    type(c_params), intent(in) :: p
    !! output directory for Gemini3D to write simulation data to (can be large files GB, TB, ...)
    integer(c_int), intent(inout) :: lid2in, lid3in  !< inout to allow optional CLI

    !> VARIABLES READ IN FROM CONFIG FILE
    real(wp) :: UTsec
    !! UT (s)
    integer, dimension(3) :: ymd
    !! year, month, day (current, not to be confused with starting year month and day in gemini_cfg structure)

    !> TEMPORAL VARIABLES
    real(wp) :: t=0._wp, dt=1e-4_wp
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
    integer :: lx1,lx2,lx3,lx2all,lx3all,lsp
    logical :: flagneuBG
    integer :: flagdneu
    real(wp) :: dtneu,dtneuBG
    integer :: myid,lid

    !> Simulation data, because these are all intended to be interoperable with C/CXX
    !    these should all be pointers, i.e. they should be allocated through specific
    !    calls and not static; this way both C and fortran main programs allocate and
    !    access these variables in analogous ways.
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    class(curvmesh), pointer :: x
    type(gemini_cfg), pointer :: cfg
    type(gemini_work), pointer :: intvars

    !> initialize message passing.  FIXME: needs to be msissetup_C()
    call mpisetup_in()
    call mpiparms(myid,lid)
    if(lid < 1) error stop 'number of MPI processes must be >= 1. Was MPI initialized properly?'

    !> command line interface
    !call cli_config_gridsize(p,lid2in,lid3in,cfg)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocations happen during this block
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cfg=>gemini_cfg_alloc()
    call cli_in(p,lid2in,lid3in,cfg)     ! transfers some data from p into cfg so cfg must be allocated prior to calling

    !> read in config file and add contents to cfg
    call read_config_in(p,cfg)           ! read configuration file and add information to cfg

    !> allocations depend on grid size so read that into our module variables
    call grid_size_in(cfg)               ! retrieve the total grid size form the input filename stored in cfg

    !> retrieve some needed module-scope variables
    call get_fullgrid_size(lx1,lx2all,lx3all)
    call get_config_vars(cfg,flagneuBG,flagdneu,dtneuBG,dtneu)

    !> MPI gridding cannot be done until we know the grid size, and needs to be done before we distribute pieces of the grid
    !    to workers
    call init_procgrid(lx2all,lx3all,lid2in,lid3in)

    !> At this point all module variables are in a state where we can set the subgrid sizes
    call calc_subgrid_size_in(lx2all,lx3all)

    !> Sizes of state variable
    call get_subgrid_size(lx1,lx2,lx3)
    call get_species_size(lsp)

    !> Allocate space for solutions, sizes will be pulled from internal modules, can happen once lx1,2,3,2all,3all defined
    !call gemini_alloc(cfg,fluidvars,fluidauxvars,electrovars,intvars)
    call gemini_double_alloc(fluidvars,fluidauxvars,electrovars)
    intvars=>gemini_work_alloc(cfg)

    !> root creates a place to put output and allocates any needed fullgrid arrays for plasma state variables
    call outdir_fullgridvaralloc(cfg,intvars,lx1,lx2all,lx3all)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> load the grid data from the input file and store in gemini module
    call read_grid_in(cfg,x)
    print*, 'Done with read_grid_in...'

    !> Set initial time variables to simulation; this requires detecting whether we are trying to restart a simulation run
    call get_initial_state(cfg,fluidvars,electrovars,intvars,x,UTsec,ymd,tdur,t,tmilestone)

    !> initialize time stepping and some aux variables
    call set_start_values_auxtimevars(it,t,tout,tglowout,tneuBG)
    call set_start_values_auxvars(x,fluidauxvars)

    !> Electric field input setup
    if(myid==0) print*, 'Priming electric field input'
    call init_Efieldinput_in(cfg,x,dt,intvars,ymd,UTsec)

    !> Recompute electrodynamic quantities needed for restarting
    !> these do not include background
    call pot2perpfield_in(x,electrovars)

    !> Get the background electric fields and compute the grid drift speed if user selected lagrangian grid, add to total field
    call BGfield_Lagrangian(cfg,x,electrovars,intvars)

    !> Precipitation input setup
    if(myid==0) print*, 'Priming precipitation input'
    call init_precipinput_in(cfg,x,dt,t,ymd,UTsec,intvars)

    !> Neutral atmosphere setup
    if(myid==0) print*, 'Computing background and priming neutral perturbation input (if used)'
    call msisinit_in(cfg)
    call init_neutralBG_in(cfg,x,dt,t,ymd,UTsec,intvars)
    call init_neutralperturb_in(dt,cfg,x,intvars,ymd,UTsec)

    !> Recompute drifts and make some decisions about whether to invoke a Lagrangian grid
    call get_initial_drifts(cfg,x,fluidvars,fluidauxvars,electrovars,intvars)

    !> control rate of console printing
    call set_update_cadence(iupdate)

    !> Main time loop
    main : do while (t < tdur)
      call dt_select(cfg,x,fluidvars,fluidauxvars,it,t,tout,tglowout,dt)

      !> get neutral background
      if ( it/=1 .and. flagneuBG .and. t>tneuBG) then              !we dont' throttle for tneuBG so we have to do things this way to not skip over...
        call cpu_time(tstart)
        call neutral_atmos_winds(cfg,x,ymd,UTsec,intvars)          ! load background states into module variables
        call neutral_atmos_wind_update(intvars)      ! apply to variables in this program unit
        tneuBG=tneuBG+dtneuBG
        if (myid==0) then
          call cpu_time(tfin)
          print *, 'Neutral background at time:  ',t,' calculated in time:  ',tfin-tstart
        end if
      end if

      !> get neutral perturbations
      if (flagdneu==1) then
        call cpu_time(tstart)
        call neutral_perturb_in(cfg,intvars,x,dt,t,ymd,UTsec)
        if (myid==0 .and. debug) then
          call cpu_time(tfin)
          print *, 'Neutral perturbations calculated in time:  ',tfin-tstart
        endif
      end if

      !> compute potential solution
      call cpu_time(tstart)
      call electrodynamics_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,it,t,dt,ymd,UTsec)
      if (myid==0 .and. debug) then
        call cpu_time(tfin)
        print *, 'Electrodynamics total solve time:  ',tfin-tstart
      endif

      !> update fluid variables
      if (myid==0 .and. debug) call cpu_time(tstart)
      call fluid_adv(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,t,dt,ymd,UTsec,(it==1),lsp,myid)
      if (myid==0 .and. debug) then
        call cpu_time(tfin)
        print *, 'Multifluid total solve time:  ',tfin-tstart
      endif

      !> Sanity check key variables before advancing
      ! FIXME: for whatever reason, it is just a fact that vs1 has trash in ghost cells after fluid_adv; I don't know why...
      call check_finite_output_in(cfg,fluidvars,electrovars,t)

      !> update time variables
      it = it + 1
      t = t + dt
      if (myid==0 .and. debug) print *, 'Moving on to time step (in sec):  ',t,'; end time of simulation:  ',tdur
      call dateinc_in(dt,ymd,UTsec)
      if (myid==0 .and. (modulo(it, iupdate) == 0 .or. debug)) then
        !! print every 10th time step to avoid extreme amounts of console printing
        print '(A,I4,A1,I0.2,A1,I0.2,A1,F12.6,A5,F8.6)', 'Current time ',ymd(1),'-',ymd(2),'-',ymd(3),' ',UTsec,'; dt=',dt
      endif

      !> see if we are doing a dry run and exit program if so
      call check_dryrun(cfg)

      !> File output
      call check_fileoutput(cfg,fluidvars,electrovars,intvars,t,tout,tglowout,tmilestone,flagoutput,ymd,UTsec)
    end do main

    !> deallocate variables and module data
    call clear_dneu_in(intvars)
    !call gemini_dealloc(cfg,fluidvars,fluidauxvars,electrovars,intvars)
    call gemini_double_dealloc(fluidvars,fluidauxvars,electrovars)
    call gemini_work_dealloc(cfg,intvars)
    call gemini_cfg_dealloc(cfg)
  end subroutine gemini_main


  !> this advances the fluid soluation by time interval dt
  subroutine fluid_adv(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,t,dt,ymd,UTsec,first,lsp,myid)
    !! J1 needed for heat conduction; E1 for momentum equation
    !! THIS SUBROUTINE ADVANCES ALL OF THE FLUID VARIABLES BY TIME STEP DT.
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: t,dt
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    logical, intent(in) :: first  !< first time step
    integer, intent(in) :: lsp
    integer, intent(in) :: myid
    real(wp) :: tstart,tfin
    real(wp) :: f107,f107a
    real(wp) :: gavg,Tninf
    integer :: isub,lsub=1   ! variables for controlling subcycling of terms

    ! pull solar indices from module type
    call get_solar_indices(cfg,f107,f107a)

    ! Prior to advection substep convert velocity and temperature to momentum and enegy density (which are local to this procedure)
    call v12rhov1_in(fluidvars,fluidauxvars)
    call T2rhoe_in(fluidvars,fluidauxvars)

    ! advection substep for all species
    call cpu_time(tstart)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Old haloing code; possibly more efficient as it only haloes one ghost cell for interface velocities
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! call halo_interface_vels_allspec_in(x,fluidvars,lsp)
    ! call interface_vels_allspec_in(fluidvars,intvars,lsp)    ! needs to happen regardless of ions v. electron due to energy eqn.
    ! call set_global_boundaries_allspec_in(x,fluidvars,fluidauxvars,intvars,lsp)
    ! call halo_allparams_in(x,fluidvars,fluidauxvars)

    ! New haloing code; probably very little performance penalty here
    call set_global_boundaries_allspec_in(x,fluidvars,fluidauxvars,intvars,lsp)
    call halo_fluidvars_in(x,fluidvars,fluidauxvars)
    call interface_vels_allspec_in(fluidvars,intvars,lsp)    ! needs to happen regardless of ions v. electron due to energy eqn.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call sweep3_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    call sweep1_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    call halo_allparams_in(x,fluidvars,fluidauxvars)
    call sweep2_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    call rhov12v1_in(fluidvars,fluidauxvars)
    call cpu_time(tfin)
    if (myid==0 .and. debug) then
      print *, 'Completed advection substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if

    ! post advection filling of null cells
    call clean_param_in(1,x,fluidvars)
    call clean_param_in(2,x,fluidvars)

    ! Compute artifical viscosity and then execute compression calculation
    call cpu_time(tstart)
    call VNRicht_artvisc_in(fluidvars,intvars)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Old haloing code; almost certainly more efficient since this only haloes one ghost cell
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! call RK2_prep_mpi_allspec_in(x,fluidvars)     ! halos velocity so we can take a divergence without artifacts

    ! This code is more general but does waste time haloing unneeded parameters and ghost cells
    call halo_fluidvars_in(x,fluidvars,fluidauxvars)
    call RK2_global_boundary_allspec_in(x,fluidvars)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call compression_in(fluidvars,fluidauxvars,intvars,x,dt)   ! this applies compression substep and then converts back to temperature
    call rhoe2T_in(fluidvars,fluidauxvars)
    call clean_param_in(3,x,fluidvars)
    call cpu_time(tfin)
    if (myid==0 .and. debug) then
      print *, 'Completed compression substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if

    ! Energy diffusion (thermal conduction) substep
    do isub=1,lsub
      ! FIXME: try to handle diffusion and sources together in call below
      !call cpu_time(tstart)
      !call energy_diffusion_in(cfg,x,fluidvars,electrovars,intvars,dt/lsub)
      !call cpu_time(tfin)
      !if (myid==0 .and. debug) then
      !  print *, 'Completed energy diffusion substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
      !end if
  
    ! cleanup and convert to specific internal energy density for sources substeps
    call clean_param_in(3,x,fluidvars)
    call T2rhoe_in(fluidvars,fluidauxvars)

    !> all workers need to "agree" on a gravity and exospheric temperature
    call get_gavg_Tinf_in(intvars,gavg,Tninf)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !> solve all source/loss processes
    call source_loss_allparams_in(cfg,fluidvars,fluidauxvars,electrovars,intvars, &
                                    x,dt/lsub,t,ymd,UTsec,f107a,f107,first,gavg,Tninf)

    ! density to be cleaned after source/loss
    call clean_param_in(3,x,fluidvars)
    call clean_param_in(2,x,fluidvars)
    call clean_param_in(1,x,fluidvars)

    end do
    !should the electron velocity be recomputed here now that densities have changed...
  end subroutine fluid_adv
end program

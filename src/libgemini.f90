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

!> This module is intended to have various interfaces/wrappers for main gemini functionality
!!   that does not involve mpi or mpi-dependent modules.  There will be a separate module with the
!!   C bindings and pointer conversions that can be called from C (i.e. wrappers for these routines).
!!   For the most part this is a bunch of "getter" routines.
module gemini3d

use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_bool, c_float, c_loc, c_null_ptr, c_ptr, c_f_pointer
use gemini_cli, only : cli
use gemini_init, only : find_config, check_input_files
use phys_consts, only: wp,debug,lnchem,lwave,lsp,pi
use meshobj, only: curvmesh
use precipdataobj, only: precipdata
use efielddataobj, only: efielddata
use neutraldataobj, only: neutraldata
use neutraldata3Dobj, only: neutraldata3D
use neutraldata3Dobj_fclaw, only: neutraldata3D_fclaw
use solfluxdataobj, only: solfluxdata
use gemini3d_config, only: gemini_cfg
use collisions, only: conductivities
use filesystem, only : expanduser
use temporal, only: cflcalc
use grid, only: grid_size,lx1,lx2,lx3,lx2all,lx3all,grid_from_extents,read_size_gridcenter, get_gridcenter, &
                  grid_internaldata_ungenerate, meshobj_alloc, meshobj_dealloc, grid_internaldata_alloc, &
                  grid_internaldata_generate, get_fullgrid_lims
use gemini3d_config, only : gemini_cfg,read_configfile
use precipBCs_mod, only: init_precipinput, precipBCs_fileinput, precipBCs
use solfluxBCs_mod, only: init_solfluxinput, solfluxBCs_fileinput, solfluxBCs
use neutral, only: neutral_info,init_neutralBG_input,neutral_atmos,neutral_winds,neutral_info_alloc,neutral_info_dealloc
use multifluid, only : sweep3_allspec_mass,sweep3_allspec_momentum,sweep3_allspec_energy, &
            sweep1_allspec_mass,sweep1_allspec_momentum,sweep1_allspec_energy, &
            sweep2_allspec_mass,sweep2_allspec_momentum,sweep2_allspec_energy, &
            VNRicht_artvisc,compression, &
            energy_diffusion,impact_ionization,solar_ionization, clean_param,rhoe2T,T2rhoe, &
            rhov12v1,v12rhov1,clean_param_after_regrid,source_loss_mass,source_loss_momentum,source_loss_energy
use advec, only: interface_vels_allspec,set_global_boundaries_allspec
use timeutils, only: dateinc
use io_nompi, only: interp_file2subgrid,plasma_output_nompi
use potential_nompi, only: set_fields_test,velocities_nompi,compute_BGEfields_nompi
use geomagnetic, only: geog2geomag,ECEFspher2ENU
use interpolation, only: interp3,interp2
use calculus, only: grad3D2,grad3D3
use sanity_check, only : check_finite_output
use potentialBCs_nompi, only: potentialBCs2D_fileinput_nompi, init_Efieldinput_nompi

implicit none (type, external)
private
public :: c_params, gemini_alloc, gemini_dealloc, init_precipinput_in, &
            set_start_values_auxtimevars, set_start_timefromcfg, set_start_values_auxvars, init_neutralBG_input_in, &
            set_update_cadence, neutral_atmos_winds, get_solar_indices, &
            v12rhov1_in, T2rhoe_in, interface_vels_allspec_in, &
            sweep3_allparams_in, sweep1_allparams_in, sweep2_allparams_in, &
            sweep3_allspec_mass_in,sweep3_allspec_momentum_in,sweep3_allspec_energy_in, &
            sweep1_allspec_mass_in,sweep1_allspec_momentum_in,sweep1_allspec_energy_in, &
            sweep2_allspec_mass_in,sweep2_allspec_momentum_in,sweep2_allspec_energy_in, &                 
            rhov12v1_in, VNRicht_artvisc_in, compression_in, rhoe2T_in, clean_param_in, &
            energy_diffusion_in, &
            source_loss_allparams_in, &
            source_loss_mass_in, source_loss_momentum_in, source_loss_energy_in, &
            clear_ionization_arrays, impact_ionization_in, solar_ionization_in, &
            dateinc_in, get_subgrid_size,get_fullgrid_size,get_config_vars, get_species_size, fluidvar_pointers, &
            fluidauxvar_pointers, electrovar_pointers, gemini_work, &
            read_fullsize_gridcenter_in, &
            gemini_work_alloc, gemini_work_dealloc, gemini_cfg_alloc, cli_in, read_config_in, gemini_cfg_dealloc, &
            grid_size_in, gemini_double_alloc, gemini_double_dealloc, gemini_grid_dealloc, &
            gemini_grid_generate, gemini_grid_generate_altnull, &
            setv2v3, v2grid, v3grid, maxcfl_in, plasma_output_nompi_in, set_global_boundaries_allspec_in, &
            get_fullgrid_lims_in,get_cfg_timevars,electrodynamics_test, precip_perturb_in, interp3_in, interp2_in, &
            check_finite_output_in, solflux_perturb_in, init_solfluxinput_in, get_it, itinc, &
            set_electrodynamics_commtype, init_efieldinput_nompi_in, efield_perturb_nompi_in

!> tracking lagrangian grid (same across all subgrids)
real(wp), protected :: v2grid,v3grid

!> internal time variables (same across all subgrids)
integer, protected :: it
real(wp), public :: tneuBG=0.0

!> type encapsulating internal arrays and parameters needed by gemini.  This is basically a catch-all for any data
!    in a gemini instance that is needed to advance the solution that must be passed into numerical procedures BUt
!    doesn't conform to simple array shapes or needs to be stored on a per-instance basis rather than globally.
type gemini_work
  !> Potential and volume emission rates
  real(wp), dimension(:,:,:), pointer :: Phiall=>null()    ! full-grid potential solution.  To store previous time step value
  real(wp), dimension(:,:,:), pointer :: iver=>null()      ! integrated volume emission rate of aurora calculated by GLOW

  !> Other variables used by the fluid solvers
  real(wp), dimension(:,:,:,:), pointer :: vs1i=>null()    ! cell interface velocities for the 1,2, and 3 directions
  real(wp), dimension(:,:,:,:), pointer :: vs2i=>null()
  real(wp), dimension(:,:,:,:), pointer :: vs3i=>null()
  real(wp), dimension(:,:,:,:), pointer :: Q=>null()       ! artificial viscosity

  !> Used to pass information about electron precipitation between procedures
  integer :: lprec=2                                                            ! number of precipitating electron populations
  real(wp), dimension(:,:,:), pointer :: W0=>null(),PhiWmWm2=>null()            ! characteristic energy and total energy flux arrays
  real(wp), dimension(:,:,:,:), pointer :: PrPrecip=>null(), Prionize=>null()   ! ionization rates from precipitation and total sources
  real(wp), dimension(:,:,:), pointer :: QePrecip=>null(), Qeionize=>null()     ! electron heating rates from precip. and total
  real(wp), dimension(:,:,:,:), pointer :: Pr=>null(),Lo=>null()                ! work arrays for tracking production/loss rates for conservation laws

  !> Use to pass information about electromagnetic boundary condtions between procedures
  integer :: flagdirich
  real(wp), dimension(:,:), pointer :: Vminx1,Vmaxx1
  real(wp), dimension(:,:), pointer :: Vminx2,Vmaxx2
  real(wp), dimension(:,:), pointer :: Vminx3,Vmaxx3
  real(wp), dimension(:,:,:), pointer :: E01,E02,E03
  real(wp), dimension(:,:), pointer :: Vminx1slab,Vmaxx1slab

  !> Used to pass solar flux data between routine
  real(wp), dimension(:,:,:,:), pointer :: Iinf

  !> Neutral information for top-level gemini program
  type(neutral_info), pointer :: atmos=>null()

  !> Inputdata objects that are needed for each subgrid
  type(precipdata), pointer :: eprecip=>null()          ! input precipitation information 
  type(efielddata), pointer :: efield=>null()           ! contains input electric field data
  class(neutraldata), pointer :: atmosperturb=>null()   ! perturbations about atmospheric background; not associated by default and may never be associated
  type(solfluxdata), pointer :: solflux=>null()         ! perturbations to solar flux, e.g., from a flare or eclipse

  real(wp), dimension(:,:,:,:), pointer :: user_output=>null()
end type gemini_work


!> type for passing C-like parameters between program units
type, bind(C) :: c_params
  !! this MUST match gemini3d.h and libgemini.f90 exactly including order
  logical(c_bool) :: fortran_nml, fortran_cli, debug, dryrun
  character(kind=c_char) :: out_dir(1000)
  !! .ini [base]
  integer(c_int) :: ymd(3)
  real(kind=c_float) :: UTsec0, tdur, dtout, activ(3), tcfl, Teinf
  !! .ini
end type c_params


contains
  !> interface subroutine from which we can read in ONLY the grid sizes
  subroutine grid_size_in(cfg)
    type(gemini_cfg), intent(in) :: cfg

    call grid_size(cfg%indatsize)
  end subroutine grid_size_in


  !> interface subroutine to handle command line inputs or otherwise setup variables that would be specified
  !    from the command line.
  subroutine cli_in(p,lid2in,lid3in,cfg)
    type(c_params), intent(in) :: p
    integer, intent(inout) :: lid2in,lid3in
    type(gemini_cfg), intent(inout) :: cfg
    character(size(p%out_dir)) :: buf
    integer :: i

    if(p%fortran_cli) then
      call cli(cfg, lid2in, lid3in, debug)
    else
      buf = "" !< ensure buf has no garbage characters

      do i = 1, len(buf)
        if (p%out_dir(i) == c_null_char) exit
        buf(i:i) = p%out_dir(i)
      enddo
      cfg%outdir = expanduser(buf)

      cfg%dryrun = p%dryrun
      debug = p%debug
    endif
  end subroutine cli_in


  !> interface layer to find and read in the config file (we assume struct has already been allocated
  subroutine read_config_in(p,cfg)
    type(c_params), intent(in) :: p
    type(gemini_cfg), intent(inout) :: cfg

    !> read the config input file, if not passed .ini info from C++ frontend
    if(p%fortran_nml) then
      call find_config(cfg)
      call read_configfile(cfg, verbose=.false.)
      call check_input_files(cfg)
    endif

    !> at this point we can check the input files and make sure we have a well-formed simulation setup
    !call check_input_files(cfg)
  end subroutine read_config_in


  !> return some data from cfg that is needed in the main program
  subroutine get_config_vars(cfg,flagneuBG,flagdneu,dtneuBG,dtneu)
    type(gemini_cfg), intent(in) :: cfg
    logical, intent(inout) :: flagneuBG
    integer, intent(inout) :: flagdneu
    real(wp), intent(inout) :: dtneuBG,dtneu

    flagneuBG=cfg%flagneuBG; flagdneu=cfg%flagdneu;
    dtneuBG=cfg%dtneuBG; dtneu=cfg%dtneu;
  end subroutine get_config_vars


  !> returns the subgrid sizes *** stored in the grid module ***
  subroutine get_subgrid_size(lx1out,lx2out,lx3out)
    integer, intent(inout) :: lx1out,lx2out,lx3out

    lx1out=lx1; lx2out=lx2; lx3out=lx3;
  end subroutine get_subgrid_size


  !> return full grid extents *** stored in the grid module ***
  subroutine get_fullgrid_size(lx1out,lx2allout,lx3allout)
    integer, intent(inout) :: lx1out,lx2allout,lx3allout

    lx1out=lx1; lx2allout=lx2all; lx3allout=lx3all;
  end subroutine get_fullgrid_size


  !> return number of species *** from phys_consts module ***
  subroutine get_species_size(lspout)
    integer, intent(inout) :: lspout

    lspout=lsp
  end subroutine get_species_size


  !> return the limits of the grid to caller
  subroutine get_fullgrid_lims_in(x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax)
    real(wp), intent(inout) :: x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax

    call get_fullgrid_lims(x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax)
  end subroutine get_fullgrid_lims_in


  !> allocate space for config struct, and return a pointer
  function gemini_cfg_alloc() result(cfg)
     type(gemini_cfg), pointer :: cfg

     allocate(cfg)
  end function gemini_cfg_alloc


  !> deallocate config struct
  subroutine gemini_cfg_dealloc(cfg)
    type(gemini_cfg), pointer, intent(inout) :: cfg

    if (associated(cfg)) then
      deallocate(cfg)
      cfg=>null()
    end if
  end subroutine gemini_cfg_dealloc


  !> allocate struct for internal variables
  function gemini_work_alloc(cfg) result(intvars)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), pointer :: intvars

    !> none of this can be done unless the size variables in the grid module are set
    if (lx1<=0 .or. lx2<=0 .or. lx3<=0 .or. lx2all<=0 .or. lx3all<=0) then
      print*,  '  Malformed size from grid module:  ',lx1,lx2,lx3,lx2all,lx3all
      error stop
    end if

    allocate(intvars)

    !> neutral variables (never need to be haloed, etc.)
    allocate(intvars%atmos)
    call neutral_info_alloc(intvars%atmos)

     !> space for integrated volume emission rates (lx2,lx3,lwave)
    if (cfg%flagglow /= 0) then
      allocate(intvars%iver(lx2,lx3,lwave))
      intvars%iver = 0
    end if

    !> allocate space for some arrays needed for fluid solves, note that these arrays are not haloed; they
    !    are computed from haloed vs1,2,3 arrays
    allocate(intvars%vs1i(1:lx1+1,1:lx2,1:lx3,1:lsp))
    allocate(intvars%vs2i(1:lx1,1:lx2+1,1:lx3,1:lsp))
    allocate(intvars%vs3i(1:lx1,1:lx2,1:lx3+1,1:lsp))
    allocate(intvars%Q(1:lx1,1:lx2,1:lx3,1:lsp))
    intvars%vs1i=0._wp
    intvars%vs2i=0._wp
    intvars%vs3i=0._wp
    intvars%Q=0._wp

    allocate(intvars%W0(1:lx2,1:lx3,1:intvars%lprec))
    allocate(intvars%PhiWmWm2,mold=intvars%W0)
    intvars%W0=1e3
    intvars%PhiWmWm2=1e-5

    ! First check that our module-scope arrays are allocated before going on to calculations.  
    ! This may need to be passed in as arguments for compatibility with trees-GEMINI
    allocate(intvars%Prprecip(1:lx1,1:lx2,1:lx3,1:lsp-1))
    intvars%Prprecip(:,:,:,:)=0.0
    allocate(intvars%Qeprecip(1:lx1,1:lx2,1:lx3))
    intvars%Qeprecip(:,:,:)=0.0
    allocate(intvars%Prionize,mold=intvars%Prprecip)
    intvars%Prionize(:,:,:,:)=0.0
    allocate(intvars%Qeionize,mold=intvars%Qeprecip)
    intvars%Qeionize(:,:,:)=0.0

    allocate(intvars%Pr(1:lx1,1:lx2,1:lx3,1:lsp))
    allocate(intvars%Lo,mold=intvars%Pr)

    allocate(intvars%Vminx1(1:lx2all,1:lx3all))
    allocate(intvars%Vmaxx1,mold=intvars%Vminx1)
    allocate(intvars%Vminx2(1:lx1,1:lx3all))
    allocate(intvars%Vmaxx2,mold=intvars%Vminx2)
    allocate(intvars%Vminx3(1:lx1,1:lx2all))
    allocate(intvars%Vmaxx3,mold=intvars%Vminx3)
    allocate(intvars%E01(1:lx1,1:lx2,1:lx3))
    allocate(intvars%E02,intvars%E03,mold=intvars%E01)
    allocate(intvars%Vminx1slab(1:lx2,1:lx3))
    allocate(intvars%Vmaxx1slab,mold=intvars%Vminx1slab)

    allocate(intvars%Iinf(1:lx1,1:lx2,1:lx3,22))   ! fix hardcoded number of wavelength bins

    allocate(intvars%eprecip)
    allocate(intvars%efield)
    ! fields of intvars%atmos are allocated in neutral:neutral_info_alloc()
    allocate(intvars%solflux)

!    ! Here the user needs to allocate any custom variables they want to pass around and/or output
!    allocate(intvars%sigP(1:lx1,1:lx2,1:lx3))
!    allocate(intvars%sigH, mold=intvars%sigP)
  end function gemini_work_alloc


  !> deallocate struct for internal variables
  subroutine gemini_work_dealloc(cfg,intvars)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), pointer, intent(inout) :: intvars

    !> neutral variables (never need to be haloed, etc.)
    !print*, 'Deallocating atmospheric state variables used in GEMINI...'
    call neutral_info_dealloc(intvars%atmos)
    deallocate(intvars%atmos)

    !> space for integrated volume emission rates (lx2,lx3,lwave)
    if (cfg%flagglow /= 0) then
      !print*, 'Deallocating glow data:  '
      deallocate(intvars%iver)
    end if

    !> allocate space for some arrays needed for fluid solves, note that these arrays are not haloed; they
    !    are computed from haloed vs1,2,3 arrays
    !print*, 'Deallocating internal variables for GEMINI...'
    deallocate(intvars%vs1i)
    deallocate(intvars%vs2i)
    deallocate(intvars%vs3i)
    deallocate(intvars%Q)

    if (associated(intvars%eprecip)) deallocate(intvars%eprecip)
    if (associated(intvars%efield)) deallocate(intvars%efield)
    !call clear_dneu(intvars%atmosperturb)    ! requies mpi so omitted here?

    deallocate(intvars%Prprecip)
    deallocate(intvars%Qeprecip)
    deallocate(intvars%Prionize)
    deallocate(intvars%Qeionize)
    if(associated(intvars%iver)) deallocate(intvars%iver)

    deallocate(intvars%Pr,intvars%Lo)

    deallocate(intvars%Vminx1,intvars%Vmaxx1)
    deallocate(intvars%Vminx2,intvars%Vmaxx2)
    deallocate(intvars%Vminx3,intvars%Vmaxx3)
    deallocate(intvars%E01,intvars%E02,intvars%E03)
    deallocate(intvars%Vminx1slab,intvars%Vmaxx1slab)

    deallocate(intvars%Iinf)

    if (associated(intvars%Phiall)) deallocate(intvars%Phiall)

    ! FIXME: why are we not deallocating intvars%eprecip and intvars%efield?

!    ! Here the user *must* deallocate their custom vars
!    deallocate(intvars%sigP,intvars%sigH)

    deallocate(intvars)
  end subroutine gemini_work_dealloc


  !> allocate space for gemini state variables, bind pointers to blocks of memory; this is primarily meant
  !    to be called from a fortran main program and simply encapsulates a set of calls to other elementary
  !    allocation procedures which could alternatively be directly called from the main GEMINI app
  subroutine gemini_alloc(cfg,fluidvars,fluidauxvars,electrovars,intvars)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars
    type(gemini_work), pointer, intent(inout) :: intvars

    !> allocate floating point arrays
    call gemini_double_alloc(fluidvars,fluidauxvars,electrovars)

    !> internal work struct
    intvars=>gemini_work_alloc(cfg)
  end subroutine gemini_alloc


  !> Fortran calls to allocate floating point arrays (should only be used from fortran)
  subroutine gemini_double_alloc(fluidvars,fluidauxvars,electrovars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars

    !> one contiguous block for overall simulation data
    allocate(fluidvars(-1:lx1+2,-1:lx2+2,-1:lx3+2,5*lsp))
    !> fluid momentum and energy density variables
    allocate(fluidauxvars(-1:lx1+2,-1:lx2+2,-1:lx3+2,2*lsp+9))
    !> electrodynamic state variables (lx1,lx2,lx3)
    allocate(electrovars(-1:lx1+2,-1:lx2+2,-1:lx3+2,7))

    !> this is a safety bit of code to make sure everything starts to zero; apparently some things are not getting initialized in some cases
    fluidvars=0._wp
    fluidauxvars=0._wp
    electrovars=0._wp
  end subroutine gemini_double_alloc


  subroutine gemini_double_dealloc(fluidvars,fluidauxvars,electrovars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars

    !> ifort generates a runtime error for this if called from C; I guess memory management needs to be done on the C-side of things
    deallocate(fluidvars)
    deallocate(fluidauxvars)
    deallocate(electrovars)
  end subroutine gemini_double_dealloc


  !> subroutine to force generate of grid internal objects (grid must already be allocated)
  subroutine gemini_grid_generate(x)
    class(curvmesh), intent(inout) :: x

    call grid_internaldata_generate(x)
  end subroutine gemini_grid_generate


  !> subroutine to force generate of grid internal objects (grid must already be allocated); user-defined null altitude
  subroutine gemini_grid_generate_altnull(x,altnull)
    class(curvmesh), intent(inout) :: x
    real(wp), intent(in) :: altnull

    call grid_internaldata_generate(x,altnull)    ! same procedure, just uses the optional argument
  end subroutine gemini_grid_generate_altnull


  !> deallocate grid data
  subroutine gemini_grid_dealloc(x,xtype,xC)
    class(curvmesh), pointer, intent(inout) :: x
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC

    call grid_internaldata_ungenerate(x)    ! this both ungenerates and also deallocates data stored in grid object
    call meshobj_dealloc(x,xtype,xC)
  end subroutine gemini_grid_dealloc


  !> force a value for the lagrangian grid drift
  subroutine setv2v3(v2gridin,v3gridin)
    real(wp), intent(in) :: v2gridin,v3gridin

    v2grid=v2gridin
    v3grid=v3gridin
  end subroutine setv2v3


  !> take a block of memory and assign pointers to various pieces representing different fluid, etc. state variables
  !!   This will be called any time a gemini library procedures needs to access individual state variables.
  subroutine fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: ns
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: vs1,vs2,vs3
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: Ts

    if (.not. associated(fluidvars)) error stop ' Attempting to bind fluid state vars to unassociated memory!'

    !> main state variables for gemini (lx1+4,lx2+4,lx3+4,lsp)
    ns=>fluidvars(:,:,:,1:lsp)
    vs1=>fluidvars(:,:,:,lsp+1:2*lsp)
    vs2=>fluidvars(:,:,:,2*lsp+1:3*lsp)
    vs3=>fluidvars(:,:,:,3*lsp+1:4*lsp)
    Ts=>fluidvars(:,:,:,4*lsp+1:5*lsp)
  end subroutine


  !> bind pointers for auxiliary fluid variables to a contiguous block of memory
  subroutine fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer, intent(inout) :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    if (.not. associated(fluidauxvars)) error stop ' Attempting to bind aux fluid state vars to unassociated memory!'

    !> pointers to aliased state variables
    rhovs1=>fluidauxvars(:,:,:,1:lsp)
    rhoes=>fluidauxvars(:,:,:,lsp+1:2*lsp)

    !> MHD-like state variables used in some calculations (lx1+4,lx2+4,lx3+4,lsp)
    rhov2=>fluidauxvars(:,:,:,2*lsp+1)
    rhov3=>fluidauxvars(:,:,:,2*lsp+2)
    B1=>fluidauxvars(:,:,:,2*lsp+3)
    B2=>fluidauxvars(:,:,:,2*lsp+4)
    B3=>fluidauxvars(:,:,:,2*lsp+5)
    v1=>fluidauxvars(:,:,:,2*lsp+6)
    v2=>fluidauxvars(:,:,:,2*lsp+7)
    v3=>fluidauxvars(:,:,:,2*lsp+8)
    rhom=>fluidauxvars(:,:,:,2*lsp+9)
  end subroutine fluidauxvar_pointers


  !> bind pointers for electomagnetic state variables to a contiguous block of memory
  subroutine electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    real(wp), dimension(:,:,:), pointer, intent(inout) :: E1,E2,E3,J1,J2,J3,Phi

    if (.not. associated(electrovars)) error stop ' Attempting to bind electro state vars to unassociated memory!'

    !> electric fields, potential, and current density
    E1=>electrovars(:,:,:,1)
    E2=>electrovars(:,:,:,2)
    E3=>electrovars(:,:,:,3)
    J1=>electrovars(:,:,:,4)
    J2=>electrovars(:,:,:,5)
    J3=>electrovars(:,:,:,6)
    Phi=>electrovars(:,:,:,7)
  end subroutine electrovar_pointers


  !> deallocate state variables include double precision data arrays; only works for all compilers from fortran main programs
  !    This is a wrapper that successively calls all deallocations and it mean to be used from a fortran app.
  subroutine gemini_dealloc(cfg,fluidvars,fluidauxvars,electrovars,intvars)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars, fluidauxvars, electrovars
    type(gemini_work), pointer, intent(inout) :: intvars

    !> ifort generates a runtime error for this if called from C; I guess memory management needs to be done on the C-side of things
    call gemini_double_dealloc(fluidvars,fluidauxvars,electrovars)

    !call gemini_dealloc_nodouble(cfg,intvars)
    call gemini_work_dealloc(cfg,intvars)
  end subroutine


  !> Basic utility to have each worker dump state variable contents to a file
  subroutine plasma_output_nompi_in(cfg,ymd,UTsec,fluidvars,electrovars,identifier,x1lims,x2lims,x3lims)
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars,electrovars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    integer, intent(in) :: identifier
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    call plasma_output_nompi(cfg%outdir,cfg%flagoutput,ymd,UTsec,ns, &
                               vs1,vs2,vs3,Ts,Phi,J1,J2,J3, &
                               identifier,x1lims,x2lims,x3lims)
  end subroutine plasma_output_nompi_in


  !> interface for pulling grid center coordinates from the input file
  subroutine read_fullsize_gridcenter_in(cfg)
    type(gemini_cfg), intent(in) :: cfg

    call read_size_gridcenter(cfg%indatsize,cfg%indatgrid)
  end subroutine read_fullsize_gridcenter_in


  !> assign initial values on some auxiliary time variables
  subroutine set_start_values_auxtimevars(t,tout,tglowout)
    real(wp), intent(inout) :: t,tout,tglowout

    !> Initialize some variables need for time stepping and output
!    it = 1; t = 0; tout = t; tglowout = t; tneuBG=t
    it = 1; tout = t; tglowout = t; tneuBG=t
  end subroutine set_start_values_auxtimevars


  ! pull relevant time variables from the cfg structure
  subroutine get_cfg_timevars(cfg,tmilestone,flagneuBG,dtneuBG,flagdneu,flagoutput)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(inout) :: tmilestone
    logical, intent(inout) :: flagneuBG
    real(wp), intent(inout) :: dtneuBG
    integer, intent(inout) :: flagdneu
    integer, intent(inout) :: flagoutput

    tmilestone=0._wp   ! make sure first output is a milestone
    flagneuBG=cfg%flagneuBG
    dtneuBG=cfg%dtneuBG
    flagdneu=cfg%flagdneu
    flagoutput=cfg%flagoutput
  end subroutine get_cfg_timevars


  !> Force time values to a specific date/time and set duration based on cfg argument fields
  subroutine set_start_timefromcfg(cfg,ymd,UTsec,tdur)
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec
    real(wp), intent(inout) :: tdur

    UTsec = cfg%UTsec0
    ymd = cfg%ymd0
    tdur = cfg%tdur
  end subroutine set_start_timefromcfg


  !> set start values for some variables not specified by the input files.
  !    some care is required here because the state variable pointers are mapped;
  !    however, note that the lbound and ubound have not been set since arrays are not passed through as dummy args
  !    with specific ubound so that we need to use intrinsic calls to make sure we fill computational cells (not ghost)
  subroutine set_start_values_auxvars(x,fluidauxvars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !ROOT/WORKERS WILL ASSUME THAT THE MAGNETIC FIELDS AND PERP FLOWS START AT ZERO
    !THIS KEEPS US FROM HAVING TO HAVE FULL-GRID ARRAYS FOR THESE STATE VARS (EXCEPT
    !FOR IN OUTPUT FNS.).  IF A SIMULATIONS IS DONE WITH INERTIAL CAPACITANCE THERE
    !WILL BE A FINITE AMOUNT OF TIME FOR THE FLOWS TO 'START UP', BUT THIS SHOULDN'T
    !BE TOO MUCH OF AN ISSUE.  WE ALSO NEED TO SET THE BACKGROUND MAGNETIC FIELD STATE
    !VARIABLE HERE TO WHATEVER IS SPECIFIED IN THE GRID STRUCTURE (THESE MUST BE CONSISTENT)
    rhov2 = 0; rhov3 = 0; v2 = 0; v3 = 0; B2 = 0; B3 = 0;
    call set_magfield(x,fluidauxvars)
  end subroutine set_start_values_auxvars


  subroutine set_magfield(x,fluidauxvars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    integer :: ix1min,ix1max,ix2min,ix2max,ix3min,ix3max
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    ix1min=lbound(B1,1)+2
    ix1max=ubound(B1,1)-2
    ix2min=lbound(B1,2)+2
    ix2max=ubound(B1,2)-2
    ix3min=lbound(B1,3)+2
    ix3max=ubound(B1,3)-2
    B1(ix1min:ix1max,ix2min:ix2max,ix3min:ix3max) = x%Bmag(1:lx1,1:lx2,1:lx3)
    !! this assumes that the grid is defined s.t. the x1 direction corresponds
    !! to the magnetic field direction (hence zero B2 and B3).
  end subroutine set_magfield


  !> Wrapper for initialization of electron precipitation data
  subroutine init_precipinput_in(cfg,x,dt,t,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call init_precipinput(dt,cfg,ymd,UTsec,x,intvars%eprecip)
  end subroutine init_precipinput_in


  !> initialize electric field input data
  subroutine init_efieldinput_nompi_in(cfg,x,dt,t,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt, t
    type(gemini_work), intent(inout) :: intvars
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call init_efieldinput_nompi(dt,cfg,ymd,UTsec,x,intvars%efield)
  end subroutine init_efieldinput_nompi_in


!  !> initialization procedure needed for MSIS 2.0
!  subroutine msisinit_in(cfg)


  subroutine init_solfluxinput_in(cfg,x,dt,t,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call init_solfluxinput(dt,cfg,ymd,UTsec,x,intvars%Iinf,intvars%solflux)
  end subroutine init_solfluxinput_in


  !> initialize the neutral background information from either MSIS
  subroutine init_neutralBG_input_in(cfg,x,dt,t,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(inout) :: x
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call init_neutralBG_input(dt,cfg,ymd,UTsec,x,v2grid,v3grid,intvars%atmos)
  end subroutine init_neutralBG_input_in


  !> set update cadence for printing out diagnostic information during simulation
  subroutine set_update_cadence(iupdate)
    integer, intent(inout) :: iupdate

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
  end subroutine set_update_cadence


  !> compute background neutral density, temperature, and wind
  subroutine neutral_atmos_winds(cfg,x,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call neutral_atmos(ymd,UTsec,x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3),x%alt(1:lx1,1:lx2,1:lx3), &
                         cfg%activ,cfg%msis_version, atmos=intvars%atmos)
    call neutral_winds(ymd, UTsec, Ap=cfg%activ(3), x=x, atmos=intvars%atmos)
  end subroutine neutral_atmos_winds


  !> get solar indices from cfg struct
  subroutine get_solar_indices(cfg,f107,f107a)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(inout) :: f107,f107a

    f107=cfg%activ(2)
    f107a=cfg%activ(1)
  end subroutine get_solar_indices


  !> convert velocity to momentum density
  subroutine v12rhov1_in(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call v12rhov1(ns,vs1,rhovs1)
  end subroutine v12rhov1_in


  !> convert temperature to specific internal energy density
  subroutine T2rhoe_in(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call T2rhoe(ns,Ts,rhoes)
  end subroutine T2rhoe_in


  !> compute interface velocities once haloing has been done
  subroutine interface_vels_allspec_in(x,fluidvars,intvars,lsp)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    type(gemini_work), intent(inout) :: intvars
    integer, intent(in) :: lsp
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call interface_vels_allspec(x,vs1,vs2,vs3,intvars%vs1i,intvars%vs2i,intvars%vs3i,lsp)    ! needs to happen regardless of ions v. electron due to energy eqn.
  end subroutine interface_vels_allspec_in


  !> enforce global boundary conditions
  subroutine set_global_boundaries_allspec_in(x,fluidvars,fluidauxvars,intvars,lsp)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    integer, intent(in) :: lsp

    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    ! fix global boundaries, as needed
    call set_global_boundaries_allspec(x%flagper,ns,rhovs1,vs1,vs2,vs3,rhoes,intvars%vs1i,lsp,x)
  end subroutine set_global_boundaries_allspec_in


  !> functions for sweeping advection
  subroutine sweep3_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !call sweep3_allparams(dt,x,intvars%vs3i,ns,rhovs1,rhoes)
    call sweep3_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
    call sweep3_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
    call sweep3_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep3_allparams_in
  subroutine sweep3_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep3_allspec_mass(dt,x,intvars%vs3i,ns)
  end subroutine sweep3_allspec_mass_in
  subroutine sweep3_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep3_allspec_momentum(dt,x,intvars%vs3i,rhovs1)
  end subroutine sweep3_allspec_momentum_in
  subroutine sweep3_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep3_allspec_energy(dt,x,intvars%vs3i,rhoes)
  end subroutine sweep3_allspec_energy_in


  subroutine sweep1_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !call sweep1_allparams(dt,x,intvars%vs1i,ns,rhovs1,rhoes)  
    call sweep1_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
    call sweep1_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
    call sweep1_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep1_allparams_in
  subroutine sweep1_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep1_allspec_mass(dt,x,intvars%vs1i,ns)
  end subroutine sweep1_allspec_mass_in
  subroutine sweep1_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep1_allspec_momentum(dt,x,intvars%vs1i,rhovs1)
  end subroutine sweep1_allspec_momentum_in
  subroutine sweep1_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep1_allspec_energy(dt,x,intvars%vs1i,rhoes)
  end subroutine sweep1_allspec_energy_in


  subroutine sweep2_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !call sweep2_allparams(dt,x,intvars%vs2i,ns,rhovs1,rhoes)   
    call sweep2_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
    call sweep2_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
    call sweep2_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep2_allparams_in
  subroutine sweep2_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep2_allspec_mass(dt,x,intvars%vs2i,ns)
  end subroutine sweep2_allspec_mass_in
  subroutine sweep2_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep2_allspec_momentum(dt,x,intvars%vs2i,rhovs1)
  end subroutine sweep2_allspec_momentum_in
  subroutine sweep2_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep2_allspec_energy(dt,x,intvars%vs2i,rhoes)
  end subroutine sweep2_allspec_energy_in


  !> conversion of momentum density to velocity
  subroutine rhov12v1_in(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call rhov12v1(ns,rhovs1,vs1)
  end subroutine rhov12v1_in


  !> compute artifical viscosity
  subroutine VNRicht_artvisc_in(fluidvars,intvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    type(gemini_work), intent(inout) :: intvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call VNRicht_artvisc(ns,vs1,intvars%Q)
  end subroutine VNRicht_artvisc_in


  !> compression substep for fluid solve
  subroutine compression_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(in) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    call compression(dt,x,vs1,vs2,vs3,intvars%Q,rhoes)   ! this applies compression substep
  end subroutine compression_in


  !> convert specific internal energy density into temperature
  subroutine rhoe2T_in(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call rhoe2T(ns,rhoes,Ts)
  end subroutine rhoe2T_in


  !> deal with null cell solutions
  subroutine clean_param_in(iparm,x,fluidvars)
    integer, intent(in) :: iparm
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: parm
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    select case (iparm)
      case (1)
        parm=>ns
      case (2)
        parm=>vs1
      case (3)
        parm=>Ts
      case default
        error stop '  libgemini:clean_params_C(); invalid parameter selected'
    end select
    call clean_param(x,iparm,parm)
  end subroutine clean_param_in


  !> diffusion of energy
  subroutine energy_diffusion_in(cfg,x,fluidvars,electrovars,intvars,dt)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(in) :: intvars
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    ! depending on the value of diffsolvetype we may want to call a completely different solver
    call energy_diffusion(dt,x,ns,Ts,J1,intvars%atmos%nn,intvars%atmos%Tn,cfg%diffsolvetype,cfg%Teinf)
  end subroutine energy_diffusion_in


  !> update the precipitation inputdata if present
  subroutine precip_perturb_in(cfg,intvars,x,dt,t,ymd,UTsec)
    real(wp), intent(in) :: t,dt
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    type(gemini_work), intent(inout) :: intvars

    if (cfg%flagprecfile==1) then
      call precipBCs_fileinput(dt,t,cfg,ymd,UTsec,x,intvars%W0,intvars%PhiWmWm2,intvars%eprecip)
    else
      !! no file input specified, so just call 'regular' function
      call precipBCs(cfg,intvars%W0,intvars%PhiWmWm2)
    end if
  end subroutine precip_perturb_in


  !> update the solar flux inputdata if present
  subroutine solflux_perturb_in(cfg,intvars,x,dt,t,ymd,UTsec)
    real(wp), intent(in) :: t,dt
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    type(gemini_work), intent(inout) :: intvars

    if (cfg%flagsolfluxfile==1) then
      call solfluxBCs_fileinput(dt,t,cfg,ymd,UTsec,x,intvars%Iinf,intvars%solflux)
    else
      !! no file input specified, so just call 'regular' function
      call solfluxBCs(cfg,x,ymd,UTsec,intvars%Iinf)
    end if
  end subroutine solflux_perturb_in
  

  !> compute boundary conditions for electric field solutions
  subroutine efield_perturb_nompi_in(cfg,intvars,x,dt,t,ymd,UTsec)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(inout) :: x     ! unit vectors could be deallocated in this procedure
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    !> assign values to the background fields, etc., irrespective of whether or not we do a potential solve
    if (cfg%flagE0file==1) then
      call potentialBCs2D_fileinput_nompi(dt,t,ymd,UTsec,cfg,x,intvars%efield,intvars%Vminx1,intvars%Vmaxx1, &
                                      intvars%Vminx2,intvars%Vmaxx2,intvars%Vminx3,intvars%Vmaxx3, &
                                      intvars%E01,intvars%E02,intvars%E03,intvars%flagdirich)
! FIXME not implemented sans mpi yet...
!    else
!      call potentialBCs2D(UTsec,cfg,x,intvars%Vminx1,intvars%Vmaxx1,intvars%Vminx2,intvars%Vmaxx2,&
!                            intvars%Vminx3,intvars%Vmaxx3, &
!                            intvars%E01,intvars%E02,intvars%E03,intvars%flagdirich)     !user needs to manually swap x2 and x3 in this function, e.g. for EIA, etc.
    end if
  end subroutine efield_perturb_nompi_in


  !> source/loss numerical solutions for all state variables; calls source/loss solutions for individual state
  !    variables, which alternatively could be called from the main program instead of this routine if one needed
  !    finer-grained control over solutions.
  subroutine source_loss_allparams_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt

    call source_loss_energy_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
    call source_loss_momentum_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
    call source_loss_mass_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
  end subroutine source_loss_allparams_in


  !> Solve for plasma mass source/losses for all species
  subroutine source_loss_mass_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)    ! not needed

    call source_loss_mass(intvars%atmos%nn,intvars%atmos%vn1,intvars%atmos%vn2,intvars%atmos%vn3, &
            intvars%atmos%Tn,ns,vs1,vs2,vs3,Ts,intvars%Pr,intvars%Lo,dt,intvars%Prionize)
  end subroutine source_loss_mass_in


  !> Momentum sources, all species
  subroutine source_loss_momentum_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    call source_loss_momentum(intvars%atmos%nn,intvars%atmos%vn1,intvars%atmos%Tn,ns,vs1,vs2,vs3,Ts,E1, &
            intvars%Q,x,intvars%Pr,intvars%Lo,dt,rhovs1)
  end subroutine source_loss_momentum_in


  !> Energy sources, all species
  subroutine source_loss_energy_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    call source_loss_energy(dt,x,cfg,ns,Ts,intvars%atmos%nn,intvars%atmos%Tn,intvars%Prionize, &
            intvars%Qeionize,intvars%atmos%vn1,intvars%atmos%vn2,intvars%atmos%vn3,vs1,vs2,vs3,rhoes, &
            intvars%Pr,intvars%Lo,intvars%Q,E2,E3)
  end subroutine source_loss_energy_in


  !> Ionization and heating rates must be re-accumulated each time so initialize to zero.  The precip arrays
  !    are not cleared here because these need to persist between time steps, e.g. glow only runs every N steps as
  !    specified by the user.  As a consequence the impact_ionization procedures need to manage initialization of 
  !    intvars%Prprecip and intvars%Qeprecip
  subroutine clear_ionization_arrays(intvars)
    type(gemini_work), intent(inout) :: intvars

    intvars%Prionize=0._wp
    intvars%Qeionize=0._wp
  end subroutine clear_ionization_arrays


  !> compute impact ionization and add results to total ionization and heating rate arrays.  Results are accumulated into
  !   intvars%Prionize and intvars%Qeprecip so these must be intialized elsewhere.  
  subroutine impact_ionization_in(cfg,fluidvars,intvars,x,dt,t,ymd, &
                                        UTsec,f107a,f107,gavg,Tninf)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: f107a,f107
    real(wp), intent(in) :: gavg,Tninf
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call impact_ionization(cfg,t,dt,x,ymd,UTsec,f107a,f107,intvars%Prprecip,intvars%Qeprecip, &
            intvars%W0,intvars%PhiWmWm2,intvars%iver,ns,Ts,intvars%atmos%nn,intvars%atmos%Tn,(get_it()==1))   ! precipiting electrons
    intvars%Prionize=intvars%Prionize+intvars%Prprecip    ! we actually need to keep a copy of the ionization by particles since GLOW not called every time step
    intvars%Qeionize=intvars%Qeionize+intvars%Qeprecip
  end subroutine impact_ionization_in


!=======
!    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
!    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
!    call source_loss_allparams(dt,t,cfg,ymd,UTsec,x,E1,E2,E3,intvars%Q,f107a,f107,intvars%atmos%nn, &
!                                     intvars%atmos%vn1,intvars%atmos%vn2,intvars%atmos%vn3, &
!                                     intvars%atmos%Tn,first,ns,rhovs1,rhoes,vs1,vs2,vs3,Ts, &
!                                     intvars%iver,gavg,Tninf,intvars%eprecip, &
!                                     cfg%diffsolvetype,cfg%Teinf,J1)
!  end subroutine source_loss_allparams_in
!>>>>>>> 969ab4efaa53ca11dd996ac250d15542314487b1


  !> Compute photoionization and *add* results to intvars%Prioinize and intvars%Qeionize
  subroutine solar_ionization_in(cfg,fluidvars,intvars,x,t,ymd,UTsec,f107a,f107,gavg,Tninf)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: f107a,f107
    real(wp), intent(in) :: gavg,Tninf
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call solar_ionization(t,x,ymd,UTsec,f107a,f107,intvars%Prionize,intvars%Qeionize,ns, &
            intvars%atmos%nn,intvars%atmos%Tn,gavg,Tninf,intvars%Iinf)     ! solar and impact ionization sources
  end subroutine solar_ionization_in


  !> Manually set the root vs. worker data collection flag in the efielddata object to user-specified value.  
  subroutine set_electrodynamics_commtype(flagrootonly, intvars)
    logical, intent(in) :: flagrootonly
    type(gemini_work), intent(inout) :: intvars

    if (associated(intvars%efield)) then
      intvars%efield%flagrootonly=flagrootonly
    else
      error stop 'Setting electro communication type without first allocating efield class...'
    end if
  end subroutine set_electrodynamics_commtype


  !> For purposes of testing we just want to set some values for the electric fields and compute drifts.  
  !    In principle this is useful for doing simulations where the potential (or background field) is
  !    specified and a potential solution is not required.  
  subroutine electrodynamics_test(cfg,x,fluidvars,fluidauxvars,electrovars,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work), intent(in) :: intvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    integer :: lx1,lx2,lx3,lsp
    real(wp), dimension(:,:,:), allocatable :: sig0,sigP,sigH,sigPgrav,sigHgrav     !FIXME: use static arrays?
    real(wp), dimension(:,:,:,:), allocatable :: muP,muH,nusn
    integer :: ix1min,ix1max,ix2min,ix2max,ix3min,ix3max

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    lx1=x%lx1; lx2=x%lx2; lx3=x%lx3; lsp=size(ns,4);

    !call set_fields_test(x,E1,E2,E3)
    allocate(sig0(lx1,lx2,lx3),sigP(lx1,lx2,lx3),sigH(lx1,lx2,lx3),sigPgrav(lx1,lx2,lx3),sigHgrav(lx1,lx2,lx3))
    allocate(muP(lx1,lx2,lx3,lsp),muH(lx1,lx2,lx3,lsp),nusn(lx1,lx2,lx3,lsp))
    call conductivities(intvars%atmos%nn,intvars%atmos%Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)

    ! pointers don't carry lbound
    ix1min=lbound(E2,1)+2
    ix1max=ubound(E2,1)-2
    ix2min=lbound(E2,2)+2
    ix2max=ubound(E2,2)-2
    ix3min=lbound(E2,3)+2
    ix3max=ubound(E2,3)-2

    if (cfg%flagE0file==1) then
      call compute_BGEfields_nompi(x,intvars%E02,intvars%E03,intvars%efield)
      intvars%E01=0._wp
      E1(ix1min:ix1max,ix2min:ix2max,ix3min:ix3max)=intvars%E01
      E2(ix1min:ix1max,ix2min:ix2max,ix3min:ix3max)=intvars%E02
      E3(ix1min:ix1max,ix2min:ix2max,ix3min:ix3max)=intvars%E03
    else 
      E1=0._wp
      E2=0._wp
      E3=0._wp
    end if

    call velocities_nompi(muP,muH,nusn,E2,E3,intvars%atmos%vn2,intvars%atmos%vn3,ns,Ts,x, &
                      cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)
    deallocate(sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
  end subroutine electrodynamics_test


  !> return the maximum cfl over the grid
  subroutine maxcfl_in(fluidvars,x,dt,maxcfl)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    class(curvmesh), pointer :: x
    real(wp), intent(in) :: dt
    real(wp), intent(inout) :: maxcfl
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call cflcalc(Ts,vs1,vs2,vs3,x%dl1i,x%dl2i,x%dl3i,dt,maxcfl)
  end subroutine


  !> Interface to access trilinear interpolation routines in gemini and interpolate MAGIC data
  subroutine interp3_in(x,y,z,q,aux,xi,yi,zi,qi,auxi)
    real(wp), dimension(:), intent(in) :: x
    real(wp), dimension(:), intent(in) :: y
    real(wp), dimension(:), intent(in) :: z
    real(wp), intent(in) :: xi,yi,zi    ! single points based on forestclaw organization
    real(wp), dimension(:), intent(inout) :: qi,auxi
    real(wp), intent(in), dimension(:,:,:,:) :: q,aux
    integer :: mx,my,mz,meqn,maux,ieqn,iaux
    real(wp), dimension(1) :: xiarr,yiarr,ziarr,qiarr,auxiarr

    ! sizes
    mx=size(q,1); my=size(q,2); mz=size(q,3);
    meqn=size(q,4); maux=size(aux,4);
    xiarr=[xi]; yiarr=[yi]; ziarr=[zi];

    ! interpolate
    do ieqn=1,meqn
      qiarr=interp3(x,y,z,q(:,:,:,ieqn),xiarr,yiarr,ziarr)     !note that the MAGIC coordinatees are permuted x,y,z
      qi(ieqn)=qiarr(1)
    end do
    do iaux=1,maux
      auxiarr=interp3(x,y,z,aux(:,:,:,iaux),xiarr,yiarr,ziarr)
      auxi(iaux)=auxiarr(1)
    end do
  end subroutine interp3_in


  !> Interface to access trilinear interpolation routines in gemini and interpolate MAGIC data
  subroutine interp2_in(x,y,q,aux,rhoi,zi,qi,auxi)
    real(wp), dimension(:), intent(in) :: x
    real(wp), dimension(:), intent(in) :: y
    real(wp), intent(in) :: rhoi,zi    ! single points based on forestclaw organization
    real(wp), dimension(:), intent(inout) :: qi,auxi
    real(wp), intent(in), dimension(:,:,:) :: q,aux
    integer :: mx,my,meqn,maux,ieqn,iaux
    real(wp), dimension(1) :: xiarr,yiarr,qiarr,auxiarr

    ! sizes
    !mx=size(q,1); my=size(q,2);
    meqn=size(q,3); maux=size(aux,3);
    xiarr=[rhoi]; yiarr=[zi];

    ! interpolate
    do ieqn=1,meqn
      qiarr=interp2(x,y,q(:,:,ieqn),xiarr,yiarr)     !note that the MAGIC coordinatees are permuted x,y,z
      qi(ieqn)=qiarr(1)
    end do
    do iaux=1,maux
      auxiarr=interp2(x,y,aux(:,:,iaux),xiarr,yiarr)
      auxi(iaux)=auxiarr(1)
    end do
  end subroutine interp2_in


  !> increment date and time arrays, this is superfluous but trying to keep outward facing function calls here.
  subroutine dateinc_in(dt,ymd,UTsec)
    real(wp), intent(in) :: dt
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec

    call dateinc(dt,ymd,UTsec)

    !print*, 'Date updated to:  ',ymd,UTsec
  end subroutine dateinc_in


  ! getter and incrementer for it variable (number of time steps taken since this start or restart.  It seems
  !   unavoidable to need almost trivial procedures for this because this must be controlled from top level 
  !   program which could be C and will need fortran procedures to modify intvar fields.  
  integer function get_it()
    get_it=it
  end function get_it
  subroutine itinc()
    it=it+1
  end subroutine itinc


  !> check main state variables for finiteness
  subroutine check_finite_output_in(cfg,fluidvars,electrovars,t)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    real(wp), intent(in) :: t

    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    call check_finite_output(cfg%outdir,t,999,vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)
  end subroutine check_finite_output_in
end module gemini3d

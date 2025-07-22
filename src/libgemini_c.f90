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

!! This module contains C/CXX wrappers for functions in libgemini.
!! These routines match those in libgemini.f90 and are
!! principally meant to convert the C pointers to various data objects into fortran pointers.
!! The grid is a class pointer (pointer to polymorphic object).
!! Other polymorphic objects (neutraldata, etc.) are kept in a static
!! derived type (intvars::gemini_work) and don't need to be passes as class pointers.

module gemini3d_C

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: iso_c_binding, only : c_int, c_bool, c_loc, c_null_ptr, c_ptr, c_f_pointer, wp => C_DOUBLE

use phys_consts, only: lnchem,lwave,lsp
use grid, only: lx1,lx2,lx3, detect_gridtype
use meshobj, only: curvmesh
use meshobj_cart, only: cartmesh
use meshobj_dipole, only: dipolemesh
use precipdataobj, only: precipdata
use efielddataobj, only: efielddata
use neutraldataobj, only: neutraldata
use gemini3d_config, only: gemini_cfg
use gemini3d, only: c_params, init_precipinput_in, &
            set_start_values_auxtimevars, set_start_values_auxvars, set_start_timefromcfg, &
            init_neutralBG_input_in, set_update_cadence, neutral_atmos_winds, get_solar_indices, &
            v12rhov1_in, T2rhoe_in, interface_vels_allspec_in, &
            sweep3_allparams_in, sweep1_allparams_in, sweep2_allparams_in, &
            sweep3_allspec_mass_in,sweep3_allspec_momentum_in,sweep3_allspec_energy_in, &
            sweep1_allspec_mass_in,sweep1_allspec_momentum_in,sweep1_allspec_energy_in, &
            sweep2_allspec_mass_in,sweep2_allspec_momentum_in,sweep2_allspec_energy_in, &
            rhov12v1_in, VNRicht_artvisc_in, compression_in, rhoe2T_in, clean_param_in, &
            energy_diffusion_in, source_loss_allparams_in, &
            source_loss_mass_in, source_loss_momentum_in, source_loss_energy_in, &
            clear_ionization_arrays, impact_ionization_in, solar_ionization_in, &
            dateinc_in, get_subgrid_size,get_fullgrid_size,get_config_vars, get_species_size, fluidvar_pointers, &
            fluidauxvar_pointers, electrovar_pointers, gemini_work, &
            read_fullsize_gridcenter_in, &
            gemini_work_alloc, gemini_work_dealloc, gemini_cfg_alloc, gemini_cfg_dealloc, grid_size_in, read_config_in, &
            cli_in, gemini_grid_generate, gemini_grid_generate_altnull, &
            gemini_grid_dealloc, setv2v3, maxcfl_in, plasma_output_nompi_in, &
            set_global_boundaries_allspec_in, get_fullgrid_lims_in, get_cfg_timevars,electrodynamics_test, &
            precip_perturb_in, interp3_in, interp2_in, check_finite_output_in, get_it, itinc, &
            set_electrodynamics_commtype, init_efieldinput_nompi_in, efield_perturb_nompi_in, &
            init_solfluxinput_in, solflux_perturb_in, source_neut_in

implicit none (type, external)

public

contains
  !> set fortran object pointer dynamic type to what is indicated in objtype.  Convert C pointer using
  !>    declared static types (c_f_pointer will not work on a polymorphic object).
  function set_gridpointer_dyntype(xtype,xC) result(x)
    type(c_ptr), intent(in) :: xC
    integer(C_INT), intent(in) :: xtype
    class(curvmesh), pointer :: x
    type(cartmesh), pointer :: xcart
    type(dipolemesh), pointer :: xdipole

    select case (xtype)
      case (1)
        call c_f_pointer(xC,xcart)
        x=>xcart
      case (2)
        call c_f_pointer(xC,xdipole)
        x=>xdipole
      case default
        write(stderr, '(a,i0)'), 'ERROR:libgemini_c:set_gridpointer_dyntype:  ' // &
          'Unable to identify object type during conversion from C to Fortran class pointer:  ',xtype
        error stop
    end select
  end function set_gridpointer_dyntype


  !> NOTE: because fortran doesn't allow you to do xcart=>x where xcart is a class extension of x the C location
  !    of the grid pointer can only be determined *at the time of creation* and cannot be arbitrarily retrieved
  !    as far as I can tell.  SO there is no inverse operation to set_gridpointer_dyntype().


  !> wrapper for command line interface
  subroutine cli_in_C(p,lid2in,lid3in,cfgC) bind(C, name='cli_in_C')
    type(c_params), intent(in) :: p
    integer(C_INT), intent(inout) :: lid2in,lid3in
    type(c_ptr), intent(inout) :: cfgC
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call cli_in(p,lid2in,lid3in,cfg)
  end subroutine cli_in_C


  !> interface for reading in the config.nml file
  subroutine read_config_in_C(p,cfgC) bind(C, name='read_config_in_C')
    type(c_params), intent(in) :: p
    type(c_ptr), intent(inout) :: cfgC
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call read_config_in(p,cfg)
  end subroutine read_config_in_C


  !> interface for reading in grid sizes into fortran module variables
  subroutine grid_size_in_C(cfgC) bind(C, name='grid_size_in_C')
    type(c_ptr), intent(in) :: cfgC
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call grid_size_in(cfg)
  end subroutine grid_size_in_C


  !> allocate a fortran struct for cfg and store the address in the C pointer cfgC
  subroutine gemini_cfg_alloc_C(cfgC) bind(C, name='gemini_cfg_alloc_C')
    type(c_ptr), intent(inout) :: cfgC
    type(gemini_cfg), pointer :: cfg

    cfg=>gemini_cfg_alloc()
    cfgC=c_loc(cfg)
  end subroutine gemini_cfg_alloc_C


  !> deallocate fortran struct connected to cfgC pointer
  subroutine gemini_cfg_dealloc_C(cfgC) bind(C, name='gemini_cfg_dealloc_C')
    type(c_ptr), intent(inout) :: cfgC
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    deallocate(cfg)
    cfg=>null()
    cfgC=c_loc(cfg)     ! send back a null pointer as a precaution
  end subroutine gemini_cfg_dealloc_C


  !> return some data from cfg that is needed in the main program
  subroutine get_config_vars_C(cfgC,flagneuBG,flagdneu,dtneuBG,dtneu) bind(C, name='get_config_vars_C')
    type(c_ptr), intent(in) :: cfgC
    logical(C_BOOL), intent(inout) :: flagneuBG
    integer(C_INT), intent(inout) :: flagdneu
    real(wp), intent(inout) :: dtneuBG,dtneu

    type(gemini_cfg), pointer :: cfg
    logical :: neuBG

    neuBG = flagneuBG

    call c_f_pointer(cfgC,cfg)
    call get_config_vars(cfg, neuBG, flagdneu,dtneuBG,dtneu)

    flagneuBG = neuBG
  end subroutine get_config_vars_C


  !> returns the subgrid sizes *** stored in the grid module ***
  subroutine get_subgrid_size_C(lx1out,lx2out,lx3out) bind(C, name='get_subgrid_size_C')
    integer(C_INT), intent(inout) :: lx1out,lx2out,lx3out

    call get_subgrid_size(lx1out,lx2out,lx3out)
  end subroutine get_subgrid_size_C


  !> return full grid extents *** stored in the grid module ***
  subroutine get_fullgrid_size_C(lx1out,lx2allout,lx3allout) bind(C, name='get_fullgrid_size_C')
    integer(C_INT), intent(inout) :: lx1out,lx2allout,lx3allout

    call get_fullgrid_size(lx1out, lx2allout, lx3allout)
  end subroutine get_fullgrid_size_C


  !> return number of species *** from phys_consts module ***
  subroutine get_species_size_C(lspout) bind(C, name='get_species_size_C')
    integer(C_INT), intent(inout) :: lspout

    call get_species_size(lspout)
  end subroutine get_species_size_C


  !> return grid limits (full grid) from module
  subroutine get_fullgrid_lims_C(x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax) bind(C,name='get_fullgrid_lims_C')
    real(wp), intent(inout) :: x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax

    call get_fullgrid_lims_in(x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax)
  end subroutine get_fullgrid_lims_C


  !> allocate space for gemini state variables, bind pointers to blocks of memory specifically internal variables
  !    we assume the C main program will itself allocate the main floating point data arrays.
  subroutine gemini_work_alloc_C(cfgC,intvarsC) bind(C, name='gemini_work_alloc_C')
    type(c_ptr), intent(in) :: cfgC
    type(c_ptr), intent(inout) :: intvarsC
    type(gemini_cfg), pointer :: cfg
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    ! allocate(intvars)
    ! call gemini_alloc_nodouble(cfg,intvars)
    intvars=>gemini_work_alloc(cfg)
    intvarsC=c_loc(intvars)
  end subroutine gemini_work_alloc_C


  !> deallocate state variables
  subroutine gemini_work_dealloc_C(cfgC,intvarsC) bind(C, name='gemini_work_dealloc_C')
    type(c_ptr), intent(in) :: cfgC
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    call c_f_pointer(intvarsC,intvars)

    !> there are issues with allocating primitives variables (doubles/ints) and then deallocating
    !    when passed back and forth with C so only deallocate the derived types
    !call gemini_dealloc_nodouble(cfg,intvars)
    call gemini_work_dealloc(cfg,intvars)
  end subroutine gemini_work_dealloc_C


  !> C wrapper for procedure to get the center location of the grid from its input file
  subroutine read_fullsize_gridcenter_C(cfgC) bind(C,name='read_fullsize_gridcenter_C')
    type(c_ptr), intent(in) :: cfgC
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call read_fullsize_gridcenter_in(cfg)
  end subroutine read_fullsize_gridcenter_C


  !> C wrapper to deallocate grid
  subroutine gemini_grid_dealloc_C(xtype,xC) bind(C, name='gemini_grid_dealloc_C')
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    class(curvmesh), pointer :: x

    !print*, 'gemini_grid_dealloc_C:  ',xtype
    x=>set_gridpointer_dyntype(xtype,xC)
    call gemini_grid_dealloc(x,xtype,xC)
  end subroutine gemini_grid_dealloc_C


  !> C wrapper to force generate of grid internal data quantities
  subroutine gemini_grid_generate_C(xtype,xC) bind(C, name='gemini_grid_generate_C')
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype,xC)
    call gemini_grid_generate(x)
  end subroutine gemini_grid_generate_C


  !> C wrapper to force generate of grid internal data quantities
  subroutine gemini_grid_generate_altnull_C(xtype,xC,altnullC) bind(C, name='gemini_grid_generate_altnull_C')
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    real(wp), intent(inout) :: altnullC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype,xC)
    call gemini_grid_generate_altnull(x,altnullC)
  end subroutine gemini_grid_generate_altnull_C


  !> wrapper to have a worker dump their state var data to a file
  subroutine plasma_output_nompi_C(cfgC,ymd,UTsec,fluidvarsC,electrovarsC, &
                                     identifier,x1lims,x2lims,x3lims) bind(C,name="plasma_output_nompi_C")
    type(c_ptr), intent(in) :: cfgC
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: electrovarsC
    integer, intent(in) :: identifier
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars

    call c_f_pointer(cfgC,cfg)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call plasma_output_nompi_in(cfg,ymd,UTsec,fluidvars,electrovars,identifier,x1lims,x2lims,x3lims)
  end subroutine plasma_output_nompi_C


  !> wrapper for forcing a particular value for the grid drift velocity
  subroutine setv2v3_C(v2gridin,v3gridin) bind(C,name='setv2v3_C')
    real(wp), intent(in) :: v2gridin,v3gridin

    call setv2v3(v2gridin,v3gridin)
  end subroutine setv2v3_C


  !> set start values for some variables.
  !    some care is required here because the state variable pointers are mapped;
  !    however, note that the lbound and ubound have not been set since arrays
  !    are not passed through as dummy args
  !    with specific ubound so that we need to use intrinsic calls to make sure we fill
  !    computational cells (not ghost)
  subroutine set_start_values_auxvars_C(xtype,xC,fluidauxvarsC) bind(C, name='set_start_values_auxvars_C')
    type(c_ptr), intent(inout) :: xC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    integer(C_INT), intent(in) :: xtype
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call set_start_values_auxvars(x,fluidauxvars)
  end subroutine set_start_values_auxvars_C


  !> initialize some auxiliary time variables used internally in gemini
  subroutine set_start_values_auxtimevars_C(t,tout,tglowout)  &
                        bind(C, name='set_start_values_auxtimevars_C')
    real(wp), intent(inout) :: t,tout,tglowout

    call set_start_values_auxtimevars(t,tout,tglowout)
  end subroutine set_start_values_auxtimevars_C


  subroutine get_cfg_timevars_C(cfgC,tmilestone,flagneuBG,dtneuBG,flagdneu,flagoutput) &
                      bind(C, name='get_cfg_timevars_C')
    type(C_PTR), intent(in) :: cfgC
    real(wp), intent(inout) :: tmilestone
    logical, intent(inout) :: flagneuBG
    real(wp), intent(inout) :: dtneuBG
    integer, intent(inout) :: flagdneu
    integer, intent(inout) :: flagoutput
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call get_cfg_timevars(cfg,tmilestone,flagneuBG,dtneuBG,flagdneu,flagoutput)
  end subroutine get_cfg_timevars_C


  !> Assign start time variables based on information in the cfg structure
  subroutine set_start_timefromcfg_C(cfgC,ymd,UTsec,tdur) bind(C, name="set_start_timefromcfg_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec
    real(wp), intent(inout) :: tdur
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call set_start_timefromcfg(cfg,ymd,UTsec,tdur)
  end subroutine set_start_timefromcfg_C


  !> Wrapper for initialization of electron precipitation data
  subroutine init_precipinput_C(cfgC,xtype,xC,dt,t,ymd,UTsec,intvarsC) bind(C, name='init_precipinput_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(intvarsC,intvars)
    call init_precipinput_in(cfg,x,dt,t,ymd,UTsec,intvars)
  end subroutine init_precipinput_C


  !> fclaw electric field input
  subroutine init_efieldinput_nompi_C(cfgC,xtype,xC,dt,t,ymd,UTsec,intvarsC) bind(C, name='init_efieldinput_nompi_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(intvarsC,intvars)
    call init_efieldinput_nompi_in(cfg,x,dt,t,ymd,UTsec,intvars)
  end subroutine init_efieldinput_nompi_C


  !> initialization procedure needed for MSIS 2.0
!  subroutine msisinit_C(cfgC) bind(C, name='msisinit_C')
!    type(c_ptr), intent(in) :: cfgC
!    type(gemini_cfg), pointer :: cfg
!
!    call c_f_pointer(cfgC,cfg)
!    call msisinit_in(cfg)
!  end subroutine msisinit_C


  !> call to initialize the neutral background information
  subroutine init_neutralBG_input_C(cfgC,xtype,xC,dt,t,ymd,UTsec,intvarsC) bind(C, name='init_neutralBG_input_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt,t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x    ! so neutral module can deallocate unit vectors once used...
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(intvarsC,intvars)
    call init_neutralBG_input_in(cfg,x,dt,t,ymd,UTsec,intvars)
  end subroutine init_neutralBG_input_C


  !> call to initialize the neutral background information
  subroutine init_solfluxinput_C(cfgC,xtype,xC,dt,t,ymd,UTsec,intvarsC) bind(C, name='init_solfluxinput_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt,t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x    ! so neutral module can deallocate unit vectors once used...
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(intvarsC,intvars)
    call init_solfluxinput_in(cfg,x,dt,t,ymd,UTsec,intvars)
  end subroutine init_solfluxinput_C 


  !> set update cadence for printing out diagnostic information during simulation
  subroutine set_update_cadence_C(iupdate) bind(C, name='set_update_cadence_C')
    integer(C_INT), intent(inout) :: iupdate

    call set_update_cadence(iupdate)
  end subroutine set_update_cadence_C


  !> compute background neutral density, temperature, and wind
  subroutine neutral_atmos_winds_C(cfgC,xtype,xC,ymd,UTsec,intvarsC) bind(C, name='neutral_atmos_winds_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(intvarsC,intvars)
    call neutral_atmos_winds(cfg,x,ymd,UTsec,intvars)
  end subroutine neutral_atmos_winds_C


  !> get solar indices from cfg struct
  subroutine get_solar_indices_C(cfgC,f107,f107a) bind(C, name='get_solar_indices_C')
    type(c_ptr), intent(in) :: cfgC
    real(wp), intent(inout) :: f107,f107a

    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call get_solar_indices(cfg,f107,f107a)
  end subroutine get_solar_indices_C


  !> convert velocity to momentum density
  subroutine v12rhov1_C(fluidvarsC,fluidauxvarsC) bind(C,name='v12rhov1_C')
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call v12rhov1_in(fluidvars,fluidauxvars)
  end subroutine v12rhov1_C


  !> convert temperature to specific internal energy density
  subroutine T2rhoe_C(fluidvarsC,fluidauxvarsC) bind(C, name='T2rhoe_C')
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call T2rhoe_in(fluidvars,fluidauxvars)
  end subroutine T2rhoe_C


  !> compute interface velocities once haloing has been done
  subroutine interface_vels_allspec_C(xtype,xC,fluidvarsC,intvarsC,lsp) bind(C, name='interface_vels_allspec_C')
    integer(C_INT), intent(in) :: xtype
    type(C_PTR), intent(in) :: xC
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: lsp

    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)
    call interface_vels_allspec_in(x,fluidvars,intvars,lsp)
  end subroutine interface_vels_allspec_C


  subroutine set_global_boundaries_allspec_C(xtype,xC, fluidvarsC,fluidauxvarsC, intvarsC, &
      lsp) bind(C, name='set_global_boundaries_allspec_C')
    integer(C_INT), intent(in) :: xtype
    type(C_PTR), intent(in) :: xC
    type(C_PTR), intent(inout) :: fluidvarsC, fluidauxvarsC
    type(C_PTR), intent(inout) :: intvarsC
    integer, intent(in) :: lsp

    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars, fluidauxvars
    type(gemini_work), pointer :: intvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)

    call set_global_boundaries_allspec_in(x, fluidvars, fluidauxvars, intvars, lsp)
  end subroutine set_global_boundaries_allspec_C


  !> functions for sweeping advection
  subroutine sweep3_allparams_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name='sweep3_allparams_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC, intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep3_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep3_allparams_C
  subroutine sweep3_allspec_mass_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="sweep3_allspec_mass_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC, intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep3_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep3_allspec_mass_C
  subroutine sweep3_allspec_momentum_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="sweep3_allspec_momentum_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC, intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep3_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep3_allspec_momentum_C
  subroutine sweep3_allspec_energy_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="sweep3_allspec_energy_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC, intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep3_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep3_allspec_energy_C


  subroutine sweep1_allparams_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name='sweep1_allparams_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep1_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep1_allparams_C
  subroutine sweep1_allspec_mass_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name='sweep1_allspec_mass_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep1_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep1_allspec_mass_C
  subroutine sweep1_allspec_momentum_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name='sweep1_allspec_momentum_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep1_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep1_allspec_momentum_C
  subroutine sweep1_allspec_energy_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name='sweep1_allspec_energy_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep1_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep1_allspec_energy_C


  subroutine sweep2_allparams_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="sweep2_allparams_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep2_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep2_allparams_C
  subroutine sweep2_allspec_mass_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="sweep2_allspec_mass_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep2_allspec_mass_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep2_allspec_mass_C
  subroutine sweep2_allspec_momentum_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="sweep2_allspec_momentum_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep2_allspec_momentum_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep2_allspec_momentum_C
  subroutine sweep2_allspec_energy_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="sweep2_allspec_energy_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep2_allspec_energy_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep2_allspec_energy_C


  !> conversion of momentum density to velocity
  subroutine rhov12v1_C(fluidvarsC, fluidauxvarsC) bind(C, name="rhov12v1_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: fluidauxvarsC

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call rhov12v1_in(fluidvars,fluidauxvars)
  end subroutine rhov12v1_C


  !> compute artifical viscosity
  subroutine VNRicht_artvisc_C(fluidvarsC,intvarsC) bind(C, name="VNRicht_artvisc_C")
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: intvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)
    call VNRicht_artvisc_in(fluidvars,intvars)
  end subroutine VNRicht_artvisc_C


  !> compression substep for fluid solve
  subroutine compression_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="compression_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call compression_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine compression_C


  !> convert specific internal energy density into temperature
  subroutine rhoe2T_C(fluidvarsC,fluidauxvarsC) bind(C, name="rhoe2T_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: fluidauxvarsC

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call rhoe2T_in(fluidvars,fluidauxvars)
  end subroutine rhoe2T_C


  !> deal with null cell solutions
  subroutine clean_param_C(iparm,xtype,xC,fluidvarsC) bind(C, name="clean_param_C")
    integer(C_INT), intent(in) :: iparm
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC

    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call clean_param_in(iparm,x,fluidvars)
  end subroutine clean_param_C


  !> diffusion of energy
  subroutine energy_diffusion_C(cfgC,xtype,xC,fluidvarsC,electrovarsC,intvarsC,dt) bind(C, name="energy_diffusion_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)
    call energy_diffusion_in(cfg,x,fluidvars,electrovars,intvars,dt)
  end subroutine energy_diffusion_C


  !> source/loss numerical solutions
  subroutine source_loss_allparams_C(cfgC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,xtype,xC,dt) &
                  bind(C, name="source_loss_allparams_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)   
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)
    call source_loss_allparams_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
  end subroutine source_loss_allparams_C


  subroutine source_loss_mass_C(cfgC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,xtype,xC,dt) bind(C, name="source_loss_mass_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x
    type(gemini_cfg), pointer :: cfg
   
    call c_f_pointer(cfgC, cfg)   
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)
    call source_loss_mass_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
  end subroutine source_loss_mass_C


  subroutine source_loss_momentum_C(cfgC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,xtype,xC,dt) &
                  bind(C, name="source_loss_momentum_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)
    call source_loss_momentum_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
  end subroutine source_loss_momentum_C


  subroutine source_loss_energy_C(cfgC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,xtype,xC,dt) &
                  bind(C, name="source_loss_energy_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)
    call source_loss_energy_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt)
  end subroutine source_loss_energy_C


  subroutine clear_ionization_arrays_C(intvarsC) bind(C, name="clear_ionization_arrays_C")
    type(c_ptr), intent(in) :: intvarsC
    type(gemini_work), pointer :: intvars

    call c_f_pointer(intvarsC,intvars)
    call clear_ionization_arrays(intvars)
  end subroutine clear_ionization_arrays_C


  subroutine source_neut_C(cfgC,fluidvarsC,intvarsC,xtype,xC) bind(C,name="source_neut_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)
    call source_neut_in(cfg,fluidvars,intvars,x)
  end subroutine source_neut_C


  subroutine impact_ionization_C(cfgC,fluidvarsC,intvarsC,xtype,xC,dt,t,ymd, &
                  UTsec,f107a,f107, & !first,
                  gavg,Tninf) &
                  bind(C, name="impact_ionization_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer(C_INT), intent(in) :: ymd(3)
    real(wp), intent(in) :: UTsec,f107,f107a
    !logical, intent(in) :: first
    real(wp), intent(in) :: gavg,Tninf

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)
    call impact_ionization_in(cfg,fluidvars,intvars,x,dt,t,ymd, &
                                        UTsec,f107a,f107,gavg,Tninf)
  end subroutine impact_ionization_C


  subroutine solar_ionization_C(cfgC,fluidvarsC,intvarsC,xtype,xC,t,ymd, &
                  UTsec,f107a,f107,gavg,Tninf) &
                  bind(C, name="solar_ionization_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: t
    integer(C_INT), intent(in) :: ymd(3)
    real(wp), intent(in) :: UTsec,f107,f107a
    real(wp), intent(in) :: gavg,Tninf

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)
    call solar_ionization_in(cfg,fluidvars,intvars,x,t,ymd, &
                                        UTsec,f107a,f107,gavg,Tninf)
  end subroutine solar_ionization_C


  !> this is a way for users to control how the efield data class reads input files (root only vs. all workers)
  subroutine set_electrodynamics_commtype_C(flagrootonlyC, intvarsC)  bind(C, name="set_electrodynamics_commtype_C")
    integer(C_INT), intent(in) :: flagrootonlyC
    type(c_ptr), intent(in) :: intvarsC
    type(gemini_work), pointer :: intvars
    logical :: flagrootonly=.true.
    
    call c_f_pointer(intvarsC,intvars)
    flagrootonly=flagrootonlyC/=0
    call set_electrodynamics_commtype(flagrootonly, intvars)
  end subroutine set_electrodynamics_commtype_C


  !> call a routine to generate test, no solve electric field information
  subroutine electrodynamics_test_C(cfgC,xtype,xC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC) &
               bind(C, name="electrodynamics_test_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)

    call electrodynamics_test(cfg,x,fluidvars,fluidauxvars,electrovars,intvars)
  end subroutine electrodynamics_test_C


  !> interface for computing cfl number
  subroutine maxcfl_C(fluidvarsC,xtype,xC,dt,maxcfl) bind(C, name="maxcfl_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt
    real(wp), intent(inout) :: maxcfl
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call maxcfl_in(fluidvars,x,dt,maxcfl)
  end subroutine maxcfl_C


  !> return the type of neutral perturbation being used
  subroutine get_neutralperturb_interptype(cfgC,interptype) bind(C,name='get_neutralperturb_interptype_C')
    type(C_PTR), intent(in) :: cfgC
    integer(C_INT), intent(inout) :: interptype
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    interptype=cfg%interptype
  end subroutine get_neutralperturb_interptype


  subroutine precip_perturb_C(cfgC, intvarsC, xtype,xC, dt,t,ymd,UTsec) bind(C, name='precip_perturb_C')
    type(C_PTR), intent(in) :: cfgC
    type(C_PTR), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(C_PTR), intent(in) :: xC
    real(wp), intent(in) :: dt,t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    type(gemini_cfg), pointer :: cfg
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)

    call precip_perturb_in(cfg,intvars,x,dt,t,ymd,UTsec)
    !call precip_perturb_in(dt,t,cfg,ymd,UTsec,x,intvars)
  end subroutine precip_perturb_C


  subroutine efield_perturb_nompi_C(cfgC, intvarsC, xtype,xC, dt,t,ymd,UTsec) bind(C, name='efield_perturb_nompi_C')
    type(C_PTR), intent(in) :: cfgC
    type(C_PTR), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(C_PTR), intent(in) :: xC
    real(wp), intent(in) :: dt,t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    type(gemini_cfg), pointer :: cfg
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)

    call efield_perturb_nompi_in(cfg,intvars,x,dt,t,ymd,UTsec)
    !call precip_perturb_in(dt,t,cfg,ymd,UTsec,x,intvars)
  end subroutine efield_perturb_nompi_C


  !> update solar fluxes stored in intvars
  subroutine solflux_perturb_C(cfgC, intvarsC, xtype,xC, dt,t,ymd,UTsec) bind(C, name='solflux_perturb_C')
    type(C_PTR), intent(in) :: cfgC
    type(C_PTR), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(C_PTR), intent(in) :: xC
    real(wp), intent(in) :: dt,t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    type(gemini_cfg), pointer :: cfg
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)

    call solflux_perturb_in(cfg,intvars,x,dt,t,ymd,UTsec)
  end subroutine solflux_perturb_C


  !> call gemini's internal interpolation code; note that input data are assumed to include ghost cells
  subroutine interp3_C(xC,yC,zC,mx,my,mz,qC,auxC,meqn,maux,xi,yi,zi,qiC,auxiC,interptype) bind(C,name='interp3_C')
    type(C_PTR), intent(in) :: xC,yC,zC
    integer(C_INT), intent(in) :: mx,my,mz
    type(C_PTR), intent(in) :: qC,auxC
    integer(C_INT), intent(in) :: meqn,maux
    real(wp), intent(in) :: xi,yi,zi
    type(C_PTR), intent(inout) :: qiC,auxiC
    integer(C_INT), intent(in) :: interptype
    real(wp), dimension(:), pointer :: x,y,z
    real(wp), dimension(:,:,:,:), pointer :: q
    real(wp), dimension(:,:,:,:), pointer :: aux
    real(wp), dimension(:), pointer :: qi
    real(wp), dimension(:), pointer :: auxi

    call c_f_pointer(xC,x,[mx+4])
    call c_f_pointer(yC,y,[my+4])
    call c_f_pointer(zC,z,[mz+4])
    call c_f_pointer(qC,q,[mx+4,my+4,mz+4,meqn])
    call c_f_pointer(auxC,aux,[mx+4,my+4,mz+4,maux])
    call c_f_pointer(qiC,qi,[meqn])
    call c_f_pointer(auxiC,auxi,[maux])

    call interp3_in(x,y,z,q,aux,xi,yi,zi,qi,auxi,interptype)
  end subroutine interp3_C


  !> call gemini's internal interpolation code
  subroutine interp2_C(xC,yC,mx,my,qC,auxC,meqn,maux,rhoi,zi,qiC,auxiC) bind(C,name='interp2_C')
    type(C_PTR), intent(in) :: xC,yC
    integer(C_INT), intent(in) :: mx,my
    type(C_PTR), intent(in) :: qC,auxC
    integer(C_INT), intent(in) :: meqn,maux
    real(wp), intent(in) :: rhoi,zi
    type(C_PTR), intent(inout) :: qiC,auxiC

    real(wp), dimension(:), pointer :: x,y
    real(wp), dimension(:,:,:), pointer :: q
    real(wp), dimension(:,:,:), pointer :: aux
    real(wp), dimension(:), pointer :: qi
    real(wp), dimension(:), pointer :: auxi

    call c_f_pointer(xC,x,[mx+4])
    call c_f_pointer(yC,y,[my+4])
    call c_f_pointer(qC,q,[mx+4,my+4,meqn])
    call c_f_pointer(auxC,aux,[mx+4,my+4,maux])
    call c_f_pointer(qiC,qi,[meqn])
    call c_f_pointer(auxiC,auxi,[maux])

    call interp2_in(x,y,q,aux,rhoi,zi,qi,auxi)
  end subroutine interp2_C


  !> increment date and time arrays, this is superfluous but trying to keep outward facing function calls here.
  subroutine dateinc_C(dt,ymd,UTsec) bind(C, name="dateinc_C")
    real(wp), intent(in) :: dt
    integer(C_INT), dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec

    call dateinc_in(dt,ymd,UTsec)
  end subroutine dateinc_C


  !> getter and incrementer for it variable (number of steps since start/restart)
  subroutine get_it_C(it) bind(C, name="get_it_C")   ! possibly not needed?
    integer(c_int), intent(inout) :: it

    it=get_it()
  end subroutine get_it_C
  subroutine itinc_C() bind(C, name="itinc_C")
    call itinc()
  end subroutine itinc_C


  subroutine check_finite_output_C(cfgC, fluidvarsC, electrovarsC, t) bind(C, name='check_finite_output_C')
    type(C_PTR), intent(in) :: cfgC
    type(C_PTR), intent(in) :: fluidvarsC, electrovarsC
    real(wp), intent(in) :: t

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars, electrovars

    call c_f_pointer(cfgC, cfg)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])

    call check_finite_output_in(cfg, fluidvars, electrovars, t)
  end subroutine check_finite_output_C
end module gemini3d_C

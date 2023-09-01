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
use gemini3d, only: c_params, init_precipinput_in, msisinit_in, &
            set_start_values_auxtimevars, set_start_values_auxvars, set_start_timefromcfg, &
            init_neutralBG_in, set_update_cadence, neutral_atmos_winds, get_solar_indices, &
            v12rhov1_in, T2rhoe_in, interface_vels_allspec_in, sweep3_allparams_in, &
            sweep1_allparams_in, sweep2_allparams_in, &
            rhov12v1_in, VNRicht_artvisc_in, compression_in, rhoe2T_in, clean_param_in, &
            energy_diffusion_in, source_loss_allparams_in, &
            dateinc_in, get_subgrid_size,get_fullgrid_size,get_config_vars, get_species_size, fluidvar_pointers, &
            fluidauxvar_pointers, electrovar_pointers, gemini_work, &
            interp_file2subgrid_in,grid_from_extents_in,read_fullsize_gridcenter_in, &
            gemini_work_alloc, gemini_work_dealloc, gemini_cfg_alloc, gemini_cfg_dealloc, grid_size_in, read_config_in, &
            cli_in, gemini_grid_generate, gemini_grid_alloc, gemini_grid_dealloc, setv2v3, maxcfl_in, plasma_output_nompi_in, &
            set_global_boundaries_allspec_in, get_fullgrid_lims_in, checkE1, get_cfg_timevars,electrodynamics_test,forceZOH_all, &
            permute_fluidvars, ipermute_fluidvars, tag4refine_location, tag4refine_vperp, clean_param_after_regrid_in, &
            get_locationsi_in,set_datainow_in, get_datainow_ptr_in, swap_statevars, interp3_in, interp2_in, &
            tag4refine_diff, tag4refine_grad, tag4coarsening_diff

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
        error stop 'unable to identify object type during conversion from C to fortran class pointer'
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


  ! FIXME: obviated needs to get rid of this here and in header file
  !> C wrapper for procedure to compute a grid object given extents and fullgrid reference point.  The class
  !    pointed to by xC must already have been allocated and assigned the correct fortran dynamic type.  
  subroutine grid_from_extents_C(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,xtype,xC) bind(C,name='grid_from_extents_C')
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer(C_INT), intent(in) :: lx1wg,lx2wg,lx3wg
    integer(C_INT), intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype,xC)
    call grid_from_extents_in(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
    ! as an extra step we need to also assign a type to the grid
    xtype=detect_gridtype(x%x1,x%x2,x%x3)
  end subroutine grid_from_extents_C


  !> C wrapper to allocate grid
  subroutine gemini_grid_alloc_C(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,xtype,xC) bind(C,name='gemini_grid_alloc_C')
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    class(curvmesh), pointer :: x

    call gemini_grid_alloc(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x,xtype,xC) 
  end subroutine gemini_grid_alloc_C


  !> C wrapper to deallocate grid
  subroutine gemini_grid_dealloc_C(xtype,xC) bind(C, name='gemini_grid_dealloc_C')
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    class(curvmesh), pointer :: x

    print*, 'gemini_grid_dealloc_C:  ',xtype
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


  !> C wrapper for procedure that reads full data from input file and interpolates to loca worker subgrid
  subroutine interp_file2subgrid_C(cfgC,xtype,xC,fluidvarsC,electrovarsC) bind(C,name='interp_file2subgrid_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: electrovarsC
    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call interp_file2subgrid_in(cfg,x,fluidvars,electrovars)
  end subroutine interp_file2subgrid_C


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
  subroutine set_start_values_auxtimevars_C(it,t,tout,tglowout,tneuBG,xtype,xC,fluidauxvarsC)  &
                        bind(C, name='set_start_values_auxtimevars_C')
    integer(C_INT), intent(inout) :: it
    real(wp), intent(inout) :: t,tout,tglowout,tneuBG
    type(c_ptr), intent(inout) :: xC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    integer(C_INT), intent(in) :: xtype
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call set_start_values_auxtimevars(it,t,tout,tglowout,tneuBG)
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


  !> initialization procedure needed for MSIS 2.0
  subroutine msisinit_C(cfgC) bind(C, name='msisinit_C')
    type(c_ptr), intent(in) :: cfgC
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call msisinit_in(cfg)
  end subroutine msisinit_C


  !> call to initialize the neutral background information
  subroutine init_neutralBG_C(cfgC,xtype,xC,dt,t,ymd,UTsec,intvarsC) bind(C, name='init_neutralBG_C')
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
    call init_neutralBG_in(cfg,x,dt,t,ymd,UTsec,intvars)
  end subroutine init_neutralBG_C


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
  subroutine interface_vels_allspec_C(fluidvarsC,intvarsC,lsp) bind(C, name='interface_vels_allspec_C')
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: lsp

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)
    call interface_vels_allspec_in(fluidvars,intvars,lsp)
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


  !> deal with null cell solutions
  subroutine clean_param_after_regrid_C(iparm,xtype,xC,fluidvarsC,intvarsC) bind(C, name="clean_param_after_regrid_C")
    integer(C_INT), intent(in) :: iparm
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: intvarsC

    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars   

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)   
    call clean_param_after_regrid_in(iparm,x,fluidvars,intvars)
  end subroutine clean_param_after_regrid_C


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
  subroutine source_loss_allparams_C(cfgC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,xtype,xC,dt,t,ymd, &
                                        UTsec,f107a,f107,first,gavg,Tninf) bind(C, name="source_loss_allparams_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt,t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: f107a,f107
    logical(C_BOOL), intent(in) :: first
    real(wp), intent(in) :: gavg,Tninf

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
    call source_loss_allparams_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt,t,ymd, &
                                        UTsec,f107a,f107,logical(first),gavg,Tninf)
  end subroutine source_loss_allparams_C


  !> echo print variable min/max for checking
  subroutine checkE1_C(fluidvarsC,fluidauxvarsC,electrovarsC,locID) bind(C, name="checkE1_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    integer, intent(in) :: locID
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call checkE1(fluidvars,fluidauxvars,electrovars,locID)
  end subroutine checkE1_C


  !> call a routine to generate test
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


  !> test to force ghost cells to zero-order hold extrapolated value
  subroutine forceZOH_all_C(fluidvarsC) bind(C,name='forceZOH_all_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call forceZOH_all(fluidvars)
  end subroutine forceZOH_all_C


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


  subroutine permute_fluidvars_C(fluidvarsC) bind(C, name='permute_fluidvars_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call permute_fluidvars(fluidvars)
  end subroutine permute_fluidvars_C


  subroutine ipermute_fluidvars_C(fluidvarsC) bind(C, name='ipermute_fluidvars_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call ipermute_fluidvars(fluidvars)
  end subroutine ipermute_fluidvars_C


  !> refinement type switchyard
  subroutine tag4refine_C(xtype,xC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,refinetype,flagrefine) &
                            bind(C, name='tag4refine_C')
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC
    integer(C_INT), intent(in) :: refinetype
    logical, intent(inout) :: flagrefine
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)

    select case (refinetype)
      case (0)    ! do a refinement based on grid location (usually used for static refines)
        call tag4refine_location(x,flagrefine)
      case (1)    ! actively refine based on perpendicular velocity magnitude (e.g. indicative of auroral or neutral forcing)
        call tag4refine_vperp(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
      case (2)
        call tag4refine_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
      case (3)
        call tag4refine_grad(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
      case default
        error stop 'unrecognized refinment criteria from user...'
    end select
  end subroutine tag4refine_C


  !> top-level routine for determining mesh coarsening
  subroutine tag4coarsening_C(xtype,xC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,flagcoarsening) &
                            bind(C, name='tag4coarsening_C')
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC
    logical, intent(inout) :: flagcoarsening
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)

    call tag4coarsening_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagcoarsening)
  end subroutine tag4coarsening_C


  !> In case we need to swap state variables
  subroutine swap_statevars_C(fluidvarsC,fluidauxvarsC) bind(C, name='swap_statevars_C')
    type(c_ptr), intent(inout) :: fluidvarsC,fluidauxvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])

    call swap_statevars(fluidvars,fluidauxvars)
  end subroutine swap_statevars_C


  !> get c pointers to magnetic coordinates of mesh sites (cell centers), note that these include ghost cells
  !  Due to the way the vtu files are created in forestclaw we likely need interface values for the cell locations
  !    As such this particular routine is likely not useful and we need separate facilities to computer and return
  !    the cell edge values.  
  subroutine get_grid_magcoords_C(xtype,xC,mlonC,mlatC,altC) bind(C, name='get_grid_magcoords_C')
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: mlonC,mlatC,altC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype, xC)

    ! for now just use geographic since already stored; also note that these do include ghost cells!
    mlonC=c_loc(x%glon)
    mlatC=c_loc(x%glat)
    altC=c_loc(x%alt)    
  end subroutine get_grid_magcoords_C


  !> force computation (if needed) of interface cell locations (geographic coordinates) and return a C pointer
  subroutine get_grid_magcoordsi_C(xtype,xC,mloniC,mlatiC,altiC) bind(C, name='get_grid_magcoordsi_C')
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: mloniC,mlatiC,altiC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype, xC)

    ! call for generation of interface locations if not already done
    call x%calc_geographici()

    ! for now just use geographic since already stored; also note that these do include ghost cells!
    mloniC=c_loc(x%gloni)
    mlatiC=c_loc(x%glati)
    altiC=c_loc(x%alti)    
  end subroutine get_grid_magcoordsi_C


  !> grab pointers and size of interpolation target locations
  subroutine get_locationsi_C(intvarsC,flagallpts,zlims,xlims,ylims, &
                                zvalsC,xvalsC,yvalsC,datavalsC,lpts,lparms) bind(C,name='get_locationsi_C')
    type(C_PTR), intent(inout) :: intvarsC
    logical, intent(in) :: flagallpts
    real(wp), dimension(2), intent(in) :: zlims,xlims,ylims   
    type(C_PTR), intent(inout) :: zvalsC,xvalsC,yvalsC,datavalsC
    integer(C_INT), intent(inout) :: lpts,lparms
    real(wp), dimension(:), pointer :: zvals,xvals,yvals
    real(wp), dimension(:,:), pointer :: datavals
    type(gemini_work), pointer :: intvars

    call c_f_pointer(intvarsC,intvars)
    call get_locationsi_in(intvars,flagallpts,zlims,xlims,ylims,zvals,xvals,yvals,datavals)
    zvalsC=c_loc(zvals)
    xvalsC=c_loc(xvals)
    yvalsC=c_loc(yvals)
    datavalsC=c_loc(datavals)
    lpts=size(datavals,1)
    lparms=size(datavals,2)
    !print*, 'Internal size chk:  ',lpts,lparms,shape(datavals)
  end subroutine get_locationsi_C


  !> retrieve (again?) the pointer to the data storage location for direct-feed
  subroutine get_datainow_ptr_C(intvarsC,datavalsC) bind(C,name='get_datainow_ptr_C')
    type(C_PTR), intent(inout) :: intvarsC
    type(C_PTR), intent(inout) :: datavalsC
    real(wp), dimension(:,:), pointer :: datavals
    type(gemini_work), pointer :: intvars  

    call c_f_pointer(intvarsC,intvars)   
    call get_datainow_ptr_in(intvars,datavals)
    datavalsC=c_loc(datavals)
  end subroutine get_datainow_ptr_C


  !> inform gemini object that the data are now in place and can be rotated/copied out
  subroutine set_datainow_C(intvarsC) bind(C,name='set_datainow_C')
    type(C_PTR), intent(inout) :: intvarsC
    type(gemini_work), pointer :: intvars

    call c_f_pointer(intvarsC,intvars)          
    call set_datainow_in(intvars)
  end subroutine set_datainow_C


  !> return the type of neutral perturbation being used
  subroutine get_neutralperturb_interptype(cfgC,interptype) bind(C,name='get_neutralperturb_interptype_C')
    type(C_PTR), intent(in) :: cfgC
    integer(C_INT), intent(inout) :: interptype
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    interptype=cfg%interptype
  end subroutine get_neutralperturb_interptype 


  !> call gemini's internal interpolation code
  subroutine interp3_C(xC,yC,zC,mx,my,mz,qC,auxC,meqn,maux,xi,yi,zi,qiC,auxiC) bind(C,name='interp3_C')
    type(C_PTR), intent(in) :: xC,yC,zC
    integer(C_INT), intent(in) :: mx,my,mz
    type(C_PTR), intent(in) :: qC,auxC
    integer(C_INT), intent(in) :: meqn,maux
    real(wp), intent(in) :: xi,yi,zi
    type(C_PTR), intent(inout) :: qiC,auxiC
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

    call interp3_in(x,y,z,q,aux,xi,yi,zi,qi,auxi)
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
end module gemini3d_C

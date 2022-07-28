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
use phys_consts, only: wp,debug,lnchem,lwave,lsp
use meshobj, only: curvmesh
use precipdataobj, only: precipdata
use efielddataobj, only: efielddata
use neutraldataobj, only: neutraldata
use gemini3d_config, only: gemini_cfg
use collisions, only: conductivities
use filesystem, only : expanduser
use temporal, only: cflcalc
use grid, only: grid_size,lx1,lx2,lx3,lx2all,lx3all,grid_from_extents,read_size_gridcenter, get_gridcenter, &
                  grid_internaldata_ungenerate, meshobj_alloc, meshobj_dealloc, grid_internaldata_alloc, &
                  grid_internaldata_generate
use gemini3d_config, only : gemini_cfg,read_configfile
use precipBCs_mod, only: init_precipinput
use msis_interface, only : msisinit
use neutral, only: neutral_info,init_neutralBG,neutral_atmos,neutral_winds,neutral_info_alloc,neutral_info_dealloc
use multifluid, only : sweep3_allparams,sweep1_allparams,sweep2_allparams,source_loss_allparams,VNRicht_artvisc,compression, &
            energy_diffusion,impact_ionization,clean_param,rhoe2T,T2rhoe, &
            rhov12v1,v12rhov1
use advec, only: interface_vels_allspec
use timeutils, only: dateinc
use io, only: interp_file2subgrid

implicit none (type, external)
private
public :: c_params, gemini_alloc, gemini_dealloc, init_precipinput_in, msisinit_in, &
            set_start_values, init_neutralBG_in, set_update_cadence, neutral_atmos_winds, get_solar_indices, &
            v12rhov1_in, T2rhoe_in, interface_vels_allspec_in, sweep3_allparams_in, &
            sweep1_allparams_in, sweep2_allparams_in, &
            rhov12v1_in, VNRicht_artvisc_in, compression_in, rhoe2T_in, clean_param_in, &
            energy_diffusion_in, source_loss_allparams_in, &
            dateinc_in, get_subgrid_size,get_fullgrid_size,get_config_vars, get_species_size, fluidvar_pointers, &
            fluidauxvar_pointers, electrovar_pointers, gemini_work, & 
            interp_file2subgrid_in,grid_from_extents_in,read_fullsize_gridcenter_in, &
            gemini_work_alloc, gemini_work_dealloc, gemini_cfg_alloc, cli_in, read_config_in, gemini_cfg_dealloc, &
            grid_size_in, gemini_double_alloc, gemini_double_dealloc, gemini_grid_alloc, gemini_grid_dealloc, &
            gemini_grid_generate, setv2v3, v2grid, v3grid, maxcfl_in


!! temp file used by MSIS 2.0
character(*), parameter :: msis2_param_file = "msis20.parm"
real(wp), protected :: v2grid,v3grid


!> type encapsulating internal arrays and parameters needed by gemini.  This is basically a catch-all for any data
!    in a gemini instance that is needed to advance the solution that must be passed into numerical procedures BUt
!    doesn't conform to simple array shapes.
type gemini_work
  real(wp), dimension(:,:,:), pointer :: Phiall=>null()    !! full-grid potential solution.  To store previous time step value
  real(wp), dimension(:,:,:), pointer :: iver    !! integrated volume emission rate of aurora calculated by GLOW

  !> Other variables used by the fluid solvers
  real(wp), dimension(:,:,:,:), pointer :: vs1i
  real(wp), dimension(:,:,:,:), pointer :: vs2i
  real(wp), dimension(:,:,:,:), pointer :: vs3i
  real(wp), dimension(:,:,:,:), pointer :: Q    ! artificial viscosity

  !> Neutral information for top-level gemini program
  type(neutral_info), pointer :: atmos

  !> Inputdata objects that are needed for each subgrid
  type(precipdata), pointer :: eprecip=>null()
  type(efielddata), pointer :: efield=>null()
  class(neutraldata), pointer :: atmosperturb=>null()   ! not associated by default and may never be associated
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
    endif 

    !> at this point we can check the input files and make sure we have a well-formed simulation setup
    call check_input_files(cfg)
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

    allocate(intvars%eprecip)
    allocate(intvars%efield)
    ! neutral stuff allocated elsewhere...
  end function gemini_work_alloc


  !> deallocate struct for internal variables
  subroutine gemini_work_dealloc(cfg,intvars)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), pointer, intent(inout) :: intvars

    !> neutral variables (never need to be haloed, etc.)
    call neutral_info_dealloc(intvars%atmos)
    deallocate(intvars%atmos)

     !> space for integrated volume emission rates (lx2,lx3,lwave)
    if (cfg%flagglow /= 0) then
      deallocate(intvars%iver)
    end if

    !> allocate space for some arrays needed for fluid solves, note that these arrays are not haloed; they
    !    are computed from haloed vs1,2,3 arrays
    deallocate(intvars%vs1i)
    deallocate(intvars%vs2i)
    deallocate(intvars%vs3i)
    deallocate(intvars%Q)

    if (associated(intvars%eprecip)) deallocate(intvars%eprecip)
    if (associated(intvars%efield)) deallocate(intvars%efield)
    !call clear_dneu(intvars%atmosperturb)    ! requies mpi so omitted here?

    if (associated(intvars%Phiall)) deallocate(intvars%Phiall)

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


  !> sequence of calls to allocate space for grid object and internal variables (analogous to gemini_work_alloc)
  subroutine gemini_grid_alloc(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x,xtype,xC)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout), pointer :: x
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    integer :: ix1,ix2,ix3
    real(wp), dimension(:), allocatable :: x1,x2,x3
    real(wp) :: glonctr,glatctr

    ! retrieve data from grid module
    call get_gridcenter(glonctr,glatctr)

    ! error checking
    if (glatctr<-90._wp .or. glatctr>90._wp) then
      error stop ' grid_from_extents:  prior to calling must use read_size_gridcenter or set_size_gridcenter to assign &
                   module variables glonctr,glatctr'
    end if

    ! create temp space
    allocate(x1(lx1wg),x2(lx2wg),x3(lx3wg))

    ! make uniformly spaced coordinate arrays
    x1=[(x1lims(1) + (x1lims(2)-x1lims(1))/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]
    x2=[(x2lims(1) + (x2lims(2)-x2lims(1))/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
    x3=[(x3lims(1) + (x3lims(2)-x3lims(1))/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)]

    ! allocate mesh class and create a C pointer to it
    call meshobj_alloc(x1,x2,x3,x,xtype,xC)

    ! allocate grid internal data arrays/structs
    call grid_internaldata_alloc(x1,x2,x3,x2,x3,glonctr,glatctr,x)

    ! get rid of temp. arrays
    deallocate(x1,x2,x3)
  end subroutine gemini_grid_alloc


  !> subroutine to force generate of grid internal objects (grid must already be allocated)
  subroutine gemini_grid_generate(x)
    class(curvmesh), intent(inout) :: x

    call grid_internaldata_generate(x)
  end subroutine gemini_grid_generate


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


  !> Have each worker read the entire input file and then interpolate it onto its own subgrid
  subroutine interp_file2subgrid_in(cfg,x,fluidvars,electrovars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars,electrovars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    call interp_file2subgrid(cfg%indatsize,cfg%indatfile,cfg%outdir,x%x1,x%x2,x%x3,ns,vs1,Ts,Phi)
  end subroutine interp_file2subgrid_in


  !> interface for pulling grid center coordinates from the input file
  subroutine read_fullsize_gridcenter_in(cfg)
    type(gemini_cfg), intent(in) :: cfg

    call read_size_gridcenter(cfg%indatsize,cfg%outdir)
  end subroutine read_fullsize_gridcenter_in


  !> A somewhat superfluous wrapper for grid generation from known extents, included here to keep with the
  !    pattern of only having applications access procedures from this module and no others.  In the case
  !    of a Cartesian grid we need to somehow set the center point so that there are glat/glon associated
  !    with each location of the mesh.  This will be obtained from the input cfg file in cases where we are
  !    generating a grid from extents.  Note that this is distinctly different from the situation where we
  !    are using read_grid() to input a grid from a file - in that case the parameters glonctr and glatctr
  !    are expected to be kept within that file.  
  subroutine grid_from_extents_in(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout) :: x

    call grid_from_extents(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)  
  end subroutine grid_from_extents_in


  !> set start values for some variables not specified by the input files.  
  !    some care is required here because the state variable pointers are mapped;
  !    however, note that the lbound and ubound have not been set since arrays are not passed through as dummy args
  !    with specific ubound so that we need to use intrinsic calls to make sure we fill computational cells (not ghost)
  subroutine set_start_values(it,t,tout,tglowout,tneuBG,x,fluidauxvars)
    integer, intent(inout) :: it
    real(wp), intent(inout) :: t,tout,tglowout,tneuBG
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    integer :: ix1min,ix1max,ix2min,ix2max,ix3min,ix3max
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !> Initialize some variables need for time stepping and output
    it = 1; t = 0; tout = t; tglowout = t; tneuBG=t

    !ROOT/WORKERS WILL ASSUME THAT THE MAGNETIC FIELDS AND PERP FLOWS START AT ZERO
    !THIS KEEPS US FROM HAVING TO HAVE FULL-GRID ARRAYS FOR THESE STATE VARS (EXCEPT
    !FOR IN OUTPUT FNS.).  IF A SIMULATIONS IS DONE WITH INERTIAL CAPACITANCE THERE
    !WILL BE A FINITE AMOUNT OF TIME FOR THE FLOWS TO 'START UP', BUT THIS SHOULDN'T
    !BE TOO MUCH OF AN ISSUE.  WE ALSO NEED TO SET THE BACKGROUND MAGNETIC FIELD STATE
    !VARIABLE HERE TO WHATEVER IS SPECIFIED IN THE GRID STRUCTURE (THESE MUST BE CONSISTENT)
    ix1min=lbound(B1,1)+2
    ix1max=ubound(B1,1)-2
    ix2min=lbound(B1,2)+2
    ix2max=ubound(B1,2)-2
    ix3min=lbound(B1,3)+2
    ix3max=ubound(B1,3)-2
    rhov2 = 0; rhov3 = 0; v2 = 0; v3 = 0; B2 = 0; B3 = 0;
    B1(ix1min:ix1max,ix2min:ix2max,ix3min:ix3max) = x%Bmag(1:lx1,1:lx2,1:lx3)
    !! this assumes that the grid is defined s.t. the x1 direction corresponds
    !! to the magnetic field direction (hence zero B2 and B3).
  end subroutine set_start_values


  !> Wrapper for initialization of electron precipitation data
  subroutine init_precipinput_in(cfg,x,dt,t,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call init_precipinput(dt,t,cfg,ymd,UTsec,x,intvars%eprecip)
  end subroutine init_precipinput_in


  !> initialization procedure needed for MSIS 2.0
  subroutine msisinit_in(cfg)
    type(gemini_cfg), intent(in) :: cfg
    logical :: exists

    if(cfg%msis_version == 20) then
      inquire(file=msis2_param_file, exist=exists)
      if(.not.exists) error stop 'could not find MSIS 2 parameter file ' // msis2_param_file // &
        ' this file must be in the same directory as gemini.bin, and run from that directory. ' // &
        'This limitation comes from how MSIS 2.x is coded internally.'
      call msisinit(parmfile=msis2_param_file)
    end if
  end subroutine msisinit_in


  !> call to initialize the neutral background information
  subroutine init_neutralBG_in(cfg,x,dt,t,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(inout) :: x    ! so neutral module can deallocate unit vectors once used...
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call init_neutralBG(dt,t,cfg,ymd,UTsec,x,v2grid,v3grid,intvars%atmos)
  end subroutine init_neutralBG_in


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

    call neutral_atmos(ymd,UTsec,x%glat,x%glon,x%alt,cfg%activ,cfg%msis_version, atmos=intvars%atmos)
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
  subroutine interface_vels_allspec_in(fluidvars,intvars,lsp)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    type(gemini_work), intent(inout) :: intvars
    integer, intent(in) :: lsp
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call interface_vels_allspec(vs1,vs2,vs3,intvars%vs1i,intvars%vs2i,intvars%vs3i,lsp)    ! needs to happen regardless of ions v. electron due to energy eqn.
  end subroutine interface_vels_allspec_in


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
    call sweep3_allparams(dt,x,intvars%vs3i,ns,rhovs1,rhoes)
  end subroutine sweep3_allparams_in

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
    call sweep1_allparams(dt,x,intvars%vs1i,ns,rhovs1,rhoes)
  end subroutine sweep1_allparams_in

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
    call sweep2_allparams(dt,x,intvars%vs2i,ns,rhovs1,rhoes)
  end subroutine sweep2_allparams_in


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
    call energy_diffusion(dt,x,ns,Ts,J1,intvars%atmos%nn,intvars%atmos%Tn,cfg%diffsolvetype,cfg%Teinf)
  end subroutine energy_diffusion_in


  !> source/loss numerical solutions
  subroutine source_loss_allparams_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt,t,ymd, &
                                        UTsec,f107a,f107,first,gavg,Tninf)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(in) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: f107a,f107
    logical, intent(in) :: first
    real(wp), intent(in) :: gavg,Tninf
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    call source_loss_allparams(dt,t,cfg,ymd,UTsec,x,E1,intvars%Q,f107a,f107,intvars%atmos%nn, &
                                     intvars%atmos%vn1,intvars%atmos%vn2,intvars%atmos%vn3, &
                                     intvars%atmos%Tn,first,ns,rhovs1,rhoes,vs1,vs2,vs3,Ts, &
                                     intvars%iver,gavg,Tninf,intvars%eprecip)
  end subroutine source_loss_allparams_in


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


  !> increment date and time arrays, this is superfluous but trying to keep outward facing function calls here.
  subroutine dateinc_in(dt,ymd,UTsec)
    real(wp), intent(in) :: dt
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec

    call dateinc(dt,ymd,UTsec)
  end subroutine dateinc_in
end module gemini3d

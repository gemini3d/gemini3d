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
!!   that does not involve mpi or mpi-dependent modules.  Note also that c bindings need to
!!   be passing only "simple" data, i.e. no structures or objects like the config struct or
!!   the grid object
module gemini3d

use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_bool, c_float, c_loc, c_null_ptr, c_ptr, c_f_pointer
use gemini_cli, only : cli
use gemini_init, only : find_config, check_input_files
use phys_consts, only: wp,debug,lnchem,lwave,lsp
use meshobj, only: curvmesh
use config, only: gemini_cfg
use collisions, only: conductivities
use filesystem, only : expanduser
use grid, only: grid_size,lx1,lx2,lx3,lx2all,lx3all
use config, only : gemini_cfg,read_configfile
use precipBCs_mod, only: init_precipinput
use msis_interface, only : msisinit
use neutral, only: init_neutralBG,neutral_atmos,neutral_winds,clear_neuBG
use multifluid, only : sweep3_allparams,sweep1_allparams,sweep2_allparams,source_loss_allparams,VNRicht_artvisc,compression, &
            energy_diffusion,impact_ionization,clean_param,rhoe2T,T2rhoe, &
            rhov12v1,v12rhov1
use advec, only: interface_vels_allspec
use timeutils, only: dateinc

implicit none (type, external)
private
public :: c_params, cli_config_gridsize, gemini_alloc, gemini_dealloc, cfg, x, init_precipinput_C, msisinit_C, &
            set_start_values, init_neutralBG_C, set_update_cadence, neutral_atmos_winds_C, get_solar_indices_C, &
            v12rhov1_C, T2rhoe_C, interface_vels_allspec_C, sweep3_allparams_C, sweep1_allparams_C, sweep2_allparams_C, &
            rhov12v1_C, VNRicht_artvisc_C, compression_C, rhoe2T_C, clean_param_C, energy_diffusion_C, source_loss_allparams_C, &
            clear_neuBG_C, dateinc_C, &
            ns,vs1,vs2,vs3,Ts,rhovs1,rhoes,E1,E2,E3,J1,J2,J3,Phi,Phiall,iver,rhov2,rhov3,B1,B2,B3,rhom,v1,v2,v3,Tn,nn,vn1, &
            vn2,vn3,vs1i,vs2i,vs3i, &
            get_subgrid_size_C,get_fullgrid_size_C,get_config_vars_C, get_species_size_C

!> these are module scope variables to avoid needing to pass as arguments in top-level main program.  In principle these could
!!   alternatively be stored in their respective modules if there is really a preference one way vs. the other.
type(gemini_cfg) :: cfg
class(curvmesh), allocatable :: x

!! Pointers to blocks of memory containing state variables
real(wp), dimension(:,:,:,:), pointer :: fluidvars
real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
real(wp), dimension(:,:,:,:), pointer :: electrovars
real(wp), dimension(:), pointer :: fluidvars_flat
real(wp), dimension(:), pointer :: fluidauxvars_flat
real(wp), dimension(:), pointer :: electrovars_flat

!! Pointers to named state variables used internally in GEMINI; these cannot be easily used in the top level program if C
!!   because mapping C-fortran pointers does not appear to work in most compilers???
real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts !! fluid state variables
real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes    ! auxiliary fluid variables for transport calculations
real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi    !! electrodynamic state variables
real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3     !! inductive state vars. (for future use - except for B1 which is used for the background field)
real(wp), dimension(:,:,:), pointer :: rhom,v1,v2,v3   !! inductive auxiliary
real(wp), dimension(:,:,:,:), pointer :: nn    !! neutral density array
real(wp), dimension(:,:,:), pointer :: Tn,vn1,vn2,vn3    !! neutral temperature and velocities
real(wp), dimension(:,:,:), pointer :: Phiall    !! full-grid potential solution.  To store previous time step value
real(wp), dimension(:,:,:), pointer :: iver    !! integrated volume emission rate of aurora calculated by GLOW

!> Other variables used by the fluid solvers
real(wp), dimension(:,:,:,:), pointer :: vs1i
real(wp), dimension(:,:,:,:), pointer :: vs2i
real(wp), dimension(:,:,:,:), pointer :: vs3i
real(wp), dimension(:,:,:,:), pointer :: Q    ! artificial viscosity

!! temp file used by MSIS 2.0
character(*), parameter :: msis2_param_file = "msis20.parm"

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
  !> basic command line and grid size determination
  subroutine cli_config_gridsize(p,lid2in,lid3in) bind(C)
    type(c_params), intent(in) :: p
    integer, intent(inout) :: lid2in,lid3in
    character(size(p%out_dir)) :: buf
    integer :: i

    !> command line interface
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

    !> read the config input file, if not passed .ini info from C++ frontend
    if(p%fortran_nml) then
      call find_config(cfg)
      call read_configfile(cfg, verbose=.false.)
    endif

    call check_input_files(cfg)

    !> read the size out of the grid file, store in module variables
    call grid_size(cfg%indatsize)
  end subroutine cli_config_gridsize


  !> return some data from cfg that is needed in the main program
  subroutine get_config_vars_C(flagneuBG,flagdneu,dtneuBG,dtneu) bind(C, name="get_config_vars_C")
    logical, intent(inout) :: flagneuBG
    integer, intent(inout) :: flagdneu
    real(wp), intent(inout) :: dtneuBG,dtneu

    flagneuBG=cfg%flagneuBG; flagdneu=cfg%flagdneu;
    dtneuBG=cfg%dtneuBG; dtneu=cfg%dtneu;
  end subroutine get_config_vars_C


  !> returns the subgrid sizes (assuming they are set to the calling procedure
  subroutine get_subgrid_size_C(lx1out,lx2out,lx3out) bind(C, name="get_subgrid_size_C")
    integer, intent(inout) :: lx1out,lx2out,lx3out

    lx1out=lx1; lx2out=lx2; lx3out=lx3;
  end subroutine


  !> return full grid extents
  subroutine get_fullgrid_size_C(lx1out,lx2allout,lx3allout) bind(C, name="get_fullgrid_size_C")
    integer, intent(inout) :: lx1out,lx2allout,lx3allout

    lx1out=lx1; lx2allout=lx2all; lx3allout=lx3all;
  end subroutine get_fullgrid_size_C


  !> return number of species from phys_consts module
  subroutine get_species_size_C(lspout) bind(C, name="get_species_size_C")
    integer, intent(inout) :: lspout

    lspout=lsp
  end subroutine get_species_size_C


  !> allocate space for gemini state variables, bind pointers to blocks of memory
  subroutine gemini_alloc(fluidvarsC,fluidauxvarsC,electrovarsC) bind(C)
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: electrovarsC

    !> one contiguous block for overall simulation data
    allocate(fluidvars(-1:lx1+2,-1:lx2+2,-1:lx3+2,5*lsp))
    call fluidvar_pointers(fluidvars)

    !> fluid momentum and energy density variables
    allocate(fluidauxvars(-1:lx1+2,-1:lx2+2,-1:lx3+2,2*lsp))
    call fluidauxvar_pointers(fluidauxvars)

    !> electrodynamic state variables (lx1,lx2,lx3)
    allocate(electrovars(-1:lx1+2,-1:lx2+2,-1:lx3+2,7))    ! include ghost cells in prep for integration with other codes
    call electrovar_pointers(electrovars)

    !> set the C pointers to the location for the memory blocks that we allocated
    fluidvarsC=c_loc(fluidvars)
    fluidauxvarsC=c_loc(fluidauxvars)
    electrovarsC=c_loc(electrovars)

    call neuMHDalloc()
  end subroutine gemini_alloc


  !> as an alternative to gemini_alloc (fortran allocation) can have C allocate space and pass in pointers that fortran will bind to state vars
  subroutine memblock_from_C(fluidvarsC,fluidauxvarsC,electrovarsC) bind(C, name="memblock_from_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: electrovarsC

    call c_f_pointer(fluidvarsC,fluidvars_flat,[(lx1+4)*(lx2+4)*(lx3+4)*(5*lsp)])
    fluidvars(-1:lx1+2,-1:lx2+2,-1:lx3+2,1:5*lsp)=>fluidvars_flat    ! this is the make absolutely sure that the bounds are okay and same as fortran alloc procedure
    call c_f_pointer(fluidauxvarsC,fluidauxvars_flat,[(lx1+4)*(lx2+4)*(lx3+4)*(2*lsp)])
    fluidauxvars(-1:lx1+2,-1:lx2+2,-1:lx3+2,1:2*lsp)=>fluidauxvars_flat
    call c_f_pointer(electrovarsC,electrovars_flat,[(lx1+4)*(lx2+4)*(lx3+4)*7])
    electrovars(-1:lx1+2,-1:lx2+2,-1:lx3+2,1:7)=>electrovars_flat

    call fluidvar_pointers(fluidvars)
    call fluidauxvar_pointers(fluidauxvars)
    call electrovar_pointers(electrovars)

    call neuMHDalloc()
  end subroutine memblock_from_C


  !> allocation space for neutral data and MHD-like parameters; this gets called regardless of whether C or Fortran allocates the main block of memory; these arrays are not visible to the "outside world"
  subroutine neuMHDalloc()
    !> MHD-like state variables used in some calculations (lx1+4,lx2+4,lx3+4,lsp)
    allocate(rhov2(-1:lx1+2,-1:lx2+2,-1:lx3+2),rhov3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
    allocate(B1(-1:lx1+2,-1:lx2+2,-1:lx3+2))
    allocate(B2(-1:lx1+2,-1:lx2+2,-1:lx3+2),B3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
    allocate(v1(-1:lx1+2,-1:lx2+2,-1:lx3+2),v2(-1:lx1+2,-1:lx2+2,-1:lx3+2), &
             v3(-1:lx1+2,-1:lx2+2,-1:lx3+2),rhom(-1:lx1+2,-1:lx2+2,-1:lx3+2))

    !> neutral variables (never need to be haloed, etc.)
    allocate(nn(lx1,lx2,lx3,lnchem),Tn(lx1,lx2,lx3),vn1(lx1,lx2,lx3), vn2(lx1,lx2,lx3),vn3(lx1,lx2,lx3))

     !> space for integrated volume emission rates (lx2,lx3,lwave)
    if (cfg%flagglow /= 0) then
      allocate(iver(lx2,lx3,lwave))
      iver = 0
    end if

    !> allocate space for some arrays needed for fluid solves, note that these arrays are not haloed; they
    !    are computed from haloed vs1,2,3 arrays.
    allocate(vs1i(1:lx1+1,1:lx2,1:lx3,1:lsp))
    allocate(vs2i(1:lx1,1:lx2+1,1:lx3,1:lsp))
    allocate(vs3i(1:lx1,1:lx2,1:lx3+1,1:lsp))
    allocate(Q(1:lx1,1:lx2,1:lx3,1:lsp))
  end subroutine 


  !> take a block of memory and assign pointers to various pieces representing different fluid, etc. state variables
  subroutine fluidvar_pointers(fluidvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars

    if (.not. associated(fluidvars)) error stop ' Attempting to bind fluid state vars to unassociated memory!'

    !> main state variables for gemini (lx1+4,lx2+4,lx3+4,lsp)
    ns=>fluidvars(:,:,:,1:lsp)
    vs1=>fluidvars(:,:,:,lsp+1:2*lsp)
    vs2=>fluidvars(:,:,:,2*lsp+1:3*lsp)
    vs3=>fluidvars(:,:,:,3*lsp+1:4*lsp)
    Ts=>fluidvars(:,:,:,4*lsp+1:5*lsp)
  end subroutine


  !> bind pointers for auxiliary fluid variables to a contiguous block of memory
  subroutine fluidauxvar_pointers(fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars

    if (.not. associated(fluidauxvars)) error stop ' Attempting to bind aux fluid state vars to unassociated memory!'

    !> pointers to aliased state variables
    rhovs1=>fluidauxvars(:,:,:,1:lsp)
    rhoes=>fluidauxvars(:,:,:,lsp+1:2*lsp)
  end subroutine fluidauxvar_pointers


  !> bind pointers for electomagnetic state variables to a contiguous block of memory
  subroutine electrovar_pointers(electrovars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars

    if (.not. associated(electrovars)) error stop ' Attempting to bind electro state vars to unassociated memory!'

    print*, 'shape of electrovars:  ', shape(electrovars)

    !> electric fields, potential, and current density
    E1=>electrovars(:,:,:,1)
    E2=>electrovars(:,:,:,2)
    E3=>electrovars(:,:,:,3)
    J1=>electrovars(:,:,:,4)
    J2=>electrovars(:,:,:,5)
    J3=>electrovars(:,:,:,6)
    Phi=>electrovars(:,:,:,7)

    print*, 'shape E2,E3:  ', shape(E2), shape(E3), minval(E2), maxval(E2), minval(E3), maxval(E3)
  end subroutine electrovar_pointers


  !> deallocate state variables
  subroutine gemini_dealloc(fluidvarsC,fluidauxvarsC,electrovarsC) bind(C)
    type(c_ptr), intent(inout) :: fluidvarsC, fluidauxvarsC, electrovarsC

    deallocate(fluidvars)
    nullify(ns,vs1,vs2,vs3,Ts)

    deallocate(fluidauxvars)
    nullify(rhovs1,rhoes)

    deallocate(rhov2,rhov3,B1,B2,B3)
    deallocate(v1,v2,v3,rhom)

    deallocate(electrovars)
    nullify(E1,E2,E3,J1,J2,J3,Phi)

    deallocate(nn,Tn,vn1,vn2,vn3)

     !> space for integrated volume emission rates
    if (cfg%flagglow /= 0) then
      deallocate(iver)
    end if

    deallocate(vs1i,vs2i,vs3i,Q)

    if (associated(Phiall)) deallocate(Phiall)

    fluidvarsC=c_null_ptr
    fluidauxvarsC=c_null_ptr
    electrovarsC=c_null_ptr
  end subroutine


  !> set start values for some variables
  subroutine set_start_values(it,t,tout,tglowout,tneuBG) bind(C)
    integer, intent(inout) :: it
    real(wp), intent(inout) :: t,tout,tglowout,tneuBG

    !> Initialize some variables need for time stepping and output
    it = 1; t = 0; tout = t; tglowout = t; tneuBG=t

    !ROOT/WORKERS WILL ASSUME THAT THE MAGNETIC FIELDS AND PERP FLOWS START AT ZERO
    !THIS KEEPS US FROM HAVING TO HAVE FULL-GRID ARRAYS FOR THESE STATE VARS (EXCEPT
    !FOR IN OUTPUT FNS.).  IF A SIMULATIONS IS DONE WITH INERTIAL CAPACITANCE THERE
    !WILL BE A FINITE AMOUNT OF TIME FOR THE FLOWS TO 'START UP', BUT THIS SHOULDN'T
    !BE TOO MUCH OF AN ISSUE.  WE ALSO NEED TO SET THE BACKGROUND MAGNETIC FIELD STATE
    !VARIABLE HERE TO WHATEVER IS SPECIFIED IN THE GRID STRUCTURE (THESE MUST BE CONSISTENT)
    rhov2 = 0; rhov3 = 0; v2 = 0; v3 = 0; B2 = 0; B3 = 0; B1(1:lx1,1:lx2,1:lx3) = x%Bmag(1:lx1,1:lx2,1:lx3)
    !! this assumes that the grid is defined s.t. the x1 direction corresponds
    !! to the magnetic field direction (hence zero B2 and B3).
  end subroutine set_start_values


  !> C binding wrapper for initialization of electron precipitation data
  subroutine init_precipinput_C(dt,t,ymd,UTsec) bind(C, name="init_precipinput_C")
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call init_precipinput(dt,t,cfg,ymd,UTsec,x)
  end subroutine init_precipinput_C


  !> initialization procedure needed for MSIS 2.0
  subroutine msisinit_C() bind(C, name="msisinit_C")
    logical :: exists

    if(cfg%msis_version == 20) then
      inquire(file=msis2_param_file, exist=exists)
      if(.not.exists) error stop 'could not find MSIS 2 parameter file ' // msis2_param_file // &
        ' this file must be in the same directory as gemini.bin, and run from that directory. ' // &
        'This limitation comes from how MSIS 2.x is coded internally.'
      call msisinit(parmfile=msis2_param_file)
    end if
  end subroutine msisinit_C


  !> call to initialize the neutral background information
  subroutine init_neutralBG_C(dt,t,ymd,UTsec,v2grid,v3grid) bind(C, name="init_neutralBG_C")
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: v2grid,v3grid

    call init_neutralBG(dt,t,cfg,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
  end subroutine init_neutralBG_C


  !> set update cadence for printing out diagnostic information during simulation
  subroutine set_update_cadence(iupdate) bind(C)
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
  subroutine neutral_atmos_winds_C(ymd,UTsec) bind(C, name="neutral_atmos_winds_C")
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call neutral_atmos(ymd,UTsec,x%glat,x%glon,x%alt,cfg%activ,cfg%msis_version)
    call neutral_winds(ymd, UTsec, Ap=cfg%activ(3), x=x)
  end subroutine neutral_atmos_winds_C


  !> get solar indices from cfg struct
  subroutine get_solar_indices_C(f107,f107a) bind(C, name="get_solar_indices_C")
    real(wp), intent(inout) :: f107,f107a

    f107=cfg%activ(2)
    f107a=cfg%activ(1)
  end subroutine get_solar_indices_C


  !> convert velocity to momentum density
  subroutine v12rhov1_C() bind(C, name="v12rhov1_C")
    call v12rhov1(ns,vs1,rhovs1)
  end subroutine v12rhov1_C


  !> convert temperature to specific internal energy density
  subroutine T2rhoe_C() bind(C, name="T2rhoe_C")
    call T2rhoe(ns,Ts,rhoes)
  end subroutine T2rhoe_C


  !> compute interface velocities once haloing has been done
  subroutine interface_vels_allspec_C(lsp) bind(C, name="interface_vels_allspec_C")
    integer, intent(in) :: lsp

    call interface_vels_allspec(vs1,vs2,vs3,vs1i,vs2i,vs3i,lsp)    ! needs to happen regardless of ions v. electron due to energy eqn.
  end subroutine interface_vels_allspec_C


  !> functions for sweeping advection
  subroutine sweep3_allparams_C(dt) bind(C, name="sweep3_allparams_C")
    real(wp), intent(in) :: dt

    call sweep3_allparams(dt,x,vs3i,ns,rhovs1,rhoes)
  end subroutine sweep3_allparams_C
  subroutine sweep1_allparams_C(dt) bind(C, name="sweep1_allparams_C")
    real(wp), intent(in) :: dt

    call sweep1_allparams(dt,x,vs1i,ns,rhovs1,rhoes)
  end subroutine sweep1_allparams_C
  subroutine sweep2_allparams_C(dt) bind(C, name="sweep2_allparams_C")
    real(wp), intent(in) :: dt

    call sweep2_allparams(dt,x,vs2i,ns,rhovs1,rhoes)
  end subroutine sweep2_allparams_C


  !> conversion of momentum density to velocity
  subroutine rhov12v1_C() bind(C, name="rhov12v1_C")
    call rhov12v1(ns,rhovs1,vs1)
  end subroutine rhov12v1_C


  !> compute artifical viscosity
  subroutine VNRicht_artvisc_C() bind(C, name="VNRicht_artvisc_C")
    call VNRicht_artvisc(ns,vs1,Q)
  end subroutine VNRicht_artvisc_C


  !> compression substep for fluid solve
  subroutine compression_C(dt) bind(C, name="compression_C")
    real(wp), intent(in) :: dt

    call compression(dt,x,vs1,vs2,vs3,Q,rhoes)   ! this applies compression substep
  end subroutine compression_C


  !> convert specific internal energy density into temperature
  subroutine rhoe2T_C() bind(C, name="rhoe2T_C")

    call rhoe2T(ns,rhoes,Ts)
  end subroutine rhoe2T_C


  !> deal with null cell solutions
  subroutine clean_param_C(iparm) bind(C, name="clean_param_C")
    integer, intent(in) :: iparm
    real(wp), dimension(:,:,:,:), pointer :: parm

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
  end subroutine clean_param_C


  !> diffusion of energy
  subroutine energy_diffusion_C(dt) bind(C, name="energy_diffusion_C")
    real(wp), intent(in) :: dt

    call energy_diffusion(dt,x,ns,Ts,J1,nn,Tn,cfg%diffsolvetype,cfg%Teinf)
  end subroutine energy_diffusion_C


  !> source/loss numerical solutions
  subroutine source_loss_allparams_C(dt,t,ymd,UTsec,f107a,f107,first,gavg,Tninf) bind(C, name="source_loss_allparams_C")
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: f107a,f107
    logical, intent(in) :: first
    real(wp), intent(in) :: gavg,Tninf

    call source_loss_allparams(dt,t,cfg,ymd,UTsec,x,E1,Q,f107a,f107,nn,vn1,vn2,vn3, &
                                     Tn,first,ns,rhovs1,rhoes,vs1,vs2,vs3,Ts,iver,gavg,Tninf)
  end subroutine source_loss_allparams_C


  !> deallocate module storage for background neutral parameters
  subroutine clear_neuBG_C() bind(C, name="clear_neuBG_C")
    call clear_neuBG()
  end subroutine clear_neuBG_C


  !> increment date and time arrays
  subroutine dateinc_C(dt,ymd,UTsec) bind(C, name="dateinc_C")
    real(wp), intent(in) :: dt
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec

    call dateinc(dt,ymd,UTsec)
  end subroutine dateinc_C
end module gemini3d

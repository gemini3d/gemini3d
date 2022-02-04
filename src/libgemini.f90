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

use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_bool, c_float
use gemini_cli, only : cli
use gemini_init, only : find_config, check_input_files
use phys_consts, only: wp,debug,lnchem,lwave,lsp
use meshobj, only: curvmesh
use config, only: gemini_cfg
use collisions, only: conductivities
use pathlib, only : expanduser
use grid, only: grid_size,lx1,lx2,lx3
use config, only : gemini_cfg,read_configfile
use precipBCs_mod, only: init_precipinput
use msis_interface, only : msisinit
use neutral, only: init_neutralBG

implicit none (type, external)
private
public :: c_params, cli_config_gridsize, gemini_alloc, gemini_dealloc, cfg, x, init_precipinput_C, msisinit_C, &
            set_start_values, init_neutralBG_C, set_update_cadence

!> these are module scope variables to avoid needing to pass as arguments in top-level main program.  In principle these could
!!   alternatively be stored in their respective modules; not sure if there is really a preference one way vs. the other.  
type(gemini_cfg) :: cfg
class(curvmesh), allocatable :: x
character(*), parameter :: msis2_param_file = "msis20.parm"

!> type for passing C-like parameters between program units
type, bind(C) :: c_params
  !! this MUST match gemini3d.h and libgemini.f90 exactly including order
  logical(c_bool) :: fortran_cli, debug, dryrun
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

    !> read the config input file 
    call find_config(cfg)
    call read_configfile(cfg, verbose=.false.)
    call check_input_files(cfg)
    
    !> read the size out of the grid file, store in module variables
    call grid_size(cfg%indatsize)
  end subroutine cli_config_gridsize


  !> allocate space for gemini state variables
  subroutine gemini_alloc(ns,vs1,vs2,vs3,Ts,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom, &
                                    E1,E2,E3,J1,J2,J3,Phi,nn,Tn,vn1,vn2,vn3,iver) bind(C)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer, intent(inout) :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom,E1,E2,E3,J1,J2,J3,Phi
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: nn
    real(wp), dimension(:,:,:), pointer, intent(inout) :: Tn,vn1,vn2,vn3
    real(wp), dimension(:,:,:), pointer, intent(inout) :: iver

    allocate(ns(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp),vs1(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp),vs2(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp), &
      vs3(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp), Ts(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp))
    allocate(rhov2(-1:lx1+2,-1:lx2+2,-1:lx3+2),rhov3(-1:lx1+2,-1:lx2+2,-1:lx3+2),B1(-1:lx1+2,-1:lx2+2,-1:lx3+2), &
             B2(-1:lx1+2,-1:lx2+2,-1:lx3+2),B3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
    allocate(v1(-1:lx1+2,-1:lx2+2,-1:lx3+2),v2(-1:lx1+2,-1:lx2+2,-1:lx3+2), &
             v3(-1:lx1+2,-1:lx2+2,-1:lx3+2),rhom(-1:lx1+2,-1:lx2+2,-1:lx3+2))
    allocate(E1(lx1,lx2,lx3),E2(lx1,lx2,lx3),E3(lx1,lx2,lx3),J1(lx1,lx2,lx3),J2(lx1,lx2,lx3),J3(lx1,lx2,lx3))
    allocate(Phi(lx1,lx2,lx3))
    allocate(nn(lx1,lx2,lx3,lnchem),Tn(lx1,lx2,lx3),vn1(lx1,lx2,lx3), vn2(lx1,lx2,lx3),vn3(lx1,lx2,lx3))

     !> space for integrated volume emission rates
    if (cfg%flagglow /= 0) then
      allocate(iver(lx2,lx3,lwave))
      iver = 0
    end if 
  end subroutine


  !> deallocate state variables
  subroutine gemini_dealloc(ns,vs1,vs2,vs3,Ts,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom,& 
                                      E1,E2,E3,J1,J2,J3,Phi,nn,Tn,vn1,vn2,vn3,iver) bind(C)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer, intent(inout) :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom,E1,E2,E3,J1,J2,J3,Phi
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: nn
    real(wp), dimension(:,:,:), pointer, intent(inout) :: Tn,vn1,vn2,vn3
    real(wp), dimension(:,:,:), pointer, intent(inout) :: iver

    deallocate(ns,vs1,vs2,vs3,Ts)
    deallocate(rhov2,rhov3,B1,B2,B3)
    deallocate(v1,v2,v3,rhom)
    deallocate(E1,E2,E3,J1,J2,J3)
    deallocate(Phi)
    deallocate(nn,Tn,vn1,vn2,vn3)

     !> space for integrated volume emission rates
    if (cfg%flagglow /= 0) then
      deallocate(iver)
    end if 
  end subroutine


  !> set start values for some variables
  subroutine set_start_values(it,t,tout,tglowout,tneuBG, &
                              rhov2,rhov3,v2,v3,B2,B3,B1) bind(C)
    integer, intent(inout) :: it
    real(wp), intent(inout) :: t,tout,tglowout,tneuBG
    real(wp), dimension(:,:,:), pointer, intent(inout) :: rhov2,rhov3,v2,v3,B2,B3,B1

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
  subroutine init_precipinput_C(dt,t,ymd,UTsec) bind(C)
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call init_precipinput(dt,t,cfg,ymd,UTsec,x)
  end subroutine init_precipinput_C


  !> initialization procedure needed for MSIS 2.0
  subroutine msisinit_C() bind(C)
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
  subroutine init_neutralBG_C(dt,t,ymd,UTsec,v2grid,v3grid,nn,Tn,vn1,vn2,vn3) bind(C)
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: v2grid,v3grid
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: nn
    real(wp), dimension(:,:,:), pointer, intent(inout) :: Tn,vn1,vn2,vn3

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
end module gemini3d

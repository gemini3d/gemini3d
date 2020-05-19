module config

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, compiler_version

use pathlib, only : expanduser, get_suffix
use phys_consts, only : wp

implicit none (type, external)
private
public :: read_configfile, gemini_cfg, get_compiler_vendor

type :: gemini_cfg

!> basic simulation information (base)
integer, dimension(3) :: ymd0
real(wp) :: UTsec0, tdur, dtout
real(wp) :: activ(3)
real(wp) :: tcfl
real(wp) :: Teinf
integer :: potsolve,flagperiodic,flagoutput
logical :: nooutput = .false.
logical :: dryrun = .false.

!> file information (files)
character(:), allocatable :: infile,outdir,indatsize,indatgrid,indatfile,out_format

!> neutral atmospheric perturbations
integer :: flagdneu                           ! whether or not to include neutral perturbations from input files
character(:), allocatable :: sourcedir        ! where the neutral input files are located
integer :: interptype=0                       ! assumptions to be used when interpolating neutrals:  0-cartesian 2D, 1-axisymmetric 2D, 3 - cartesian 3D
real(wp) :: sourcemlat=0,sourcemlon=0         ! source latitude and longitude
real(wp) :: dtneu=0                           ! time step between neutral inputs
real(wp) :: dxn=0,drhon=0,dzn=0               ! dx,dy (or drho),dz for neutral inputs

!> preciptiation file inputs
integer :: flagprecfile                       ! whether or not we have precipitation input from a file
character(:), allocatable :: precdir          ! location of precipitation input data
real(wp) :: dtprec=0                          ! time step between precipitation inputs

!> electric field file inputs
integer :: flagE0file                         ! whether or not to have electric field file input
character(:), allocatable :: E0dir            ! location of electric field input data
real(wp) :: dtE0=0                            ! time step between electric field inputs

!> GLOW parameters
integer :: flagglow              ! whether or not to use glow to compute impact ionization
real(wp) :: dtglow, dtglowout    ! time step between GLOW updates and outputs for GLOW emissions

!! parameters below this line can only be changed via the .nml input format
!> equatorial ionization anomaly
logical :: flagEIA=.false.       ! whether or not to include and equatorial ionization anomaly in simulation
real(wp) :: v0equator=10._wp      ! max vertical drift of plasma at equator for EIA

!> varying neutral atmosphere background
logical :: flagneuBG=.false.                ! whether or not to allow MSIS to be called to update neutral background
real(wp) :: dtneuBG=900._wp                  ! approximate time between MSIS calls

!> background preciptation
real(wp) :: PhiWBG=1e-3_wp                      ! background total energy flux in mW/m^2
real(wp) :: W0BG=3e3_wp                      ! background characteristic energy for precipitation

!> parallel current calculations
logical :: flagJpar=.true.                  ! whether or not to compute parallel current (some simulation setups will give really poor results); code ignores this if potential is resolved along the field line since computing Jpar will not be prone to artifacts as it is in th EFL cases...

!> inertial capacitance
integer :: flagcap = 0           ! use inertial capacitance? 0 - set all to zero, 1 - use ionosphere to compute, 2 - add a magnetospheric part
real(wp) :: magcap=30._wp        ! value of integrated magnetospheric capacitance to use

!> type of diffusion solver to sue
integer :: diffsolvetype=2       ! 1 - first order backward Euler time stepping; 2 - 2nd order TRBDF2 diffusion solver

end type gemini_cfg


interface
module subroutine read_nml(cfg, verbose)
class(gemini_cfg), intent(inout) :: cfg
logical, intent(in), optional :: verbose
end subroutine read_nml

module subroutine read_ini(cfg)
class(gemini_cfg), intent(inout) :: cfg
end subroutine read_ini
end interface


contains


subroutine read_configfile(cfg, verbose)

class(gemini_cfg), intent(inout) :: cfg
logical, intent(in), optional :: verbose

!! READS THE INPUT CONFIGURAITON FILE, ASSIGNS VARIABLES FOR FILENAMES, SIZES, ETC.
integer :: i, realbits, lxp, lyp
real(wp) :: NaN, glat, glon, xdist, ydist, alt_min, alt_max, alt_scale(4), Bincl, nmf, nme

!> READ CONFIG FILE FOR THIS SIMULATION
!! NOTE: Namelist file groups must be read in order they appear in the Namelist file, or End of File error occurs

if (cfg%infile(len(cfg%infile)-3 : len(cfg%infile)) == '.nml') then
  call read_nml(cfg, verbose)
else
  call read_ini(cfg)
endif

end subroutine read_configfile


function get_compiler_vendor() result(vendor)
character(:), allocatable :: vendor
character(80) :: cvers
integer :: i, j
character(*), parameter :: vendors(2) = [character(5) :: "Intel", "GCC"]

cvers = compiler_version()

do j = 1,size(vendors)
  i = index(cvers, trim(vendors(j)))
  if (i > 0) then
    vendor = trim(vendors(j))
    exit
  endif
enddo

if(vendor=="GCC") vendor = "GNU"

if(allocated(vendor)) return

vendor = ""
write(stderr,'(A,/,A)') "could not determine compiler vendor from",cvers

end function get_compiler_vendor

end module config

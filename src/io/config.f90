module config

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, compiler_version

use pathlib, only : expanduser, get_suffix, make_absolute
use phys_consts, only : wp

implicit none (type, external)
private
public :: read_configfile, gemini_cfg, get_compiler_vendor, expand_envvar

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
character(:), allocatable :: infile,outdir,indatsize,indatgrid,indatfile,out_format, &
  fieldpointfile

character(:), allocatable :: git_revision

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

integer :: flag_fang=2008   !< configure Fang ionization model

!! parameters below this line can only be changed via the .nml input format
!> equatorial ionization anomaly
logical :: flagEIA=.false.       ! whether or not to include and equatorial ionization anomaly in simulation
real(wp) :: v0equator=10._wp      ! max vertical drift of plasma at equator for EIA

!> varying neutral atmosphere background
logical :: flagneuBG=.false.                ! whether or not to allow MSIS to be called to update neutral background
real(wp) :: dtneuBG=900._wp                  ! approximate time between MSIS calls
integer :: msis_version

!> background preciptation
real(wp) :: PhiWBG=1e-3_wp                      ! background total energy flux in mW/m^2
real(wp) :: W0BG=3e3_wp                      ! background characteristic energy for precipitation

!> parallel current calculations
logical :: flagJpar=.true.                  ! whether or not to compute parallel current (some simulation setups will give really poor results); code ignores this if potential is resolved along the field line since computing Jpar will not be prone to artifacts as it is in th EFL cases...

!> inertial capacitance
integer :: flagcap = 0           ! use inertial capacitance? 0 - set all to zero, 1 - use ionosphere to compute, 2 - add a magnetospheric part
real(wp) :: magcap=5._wp        ! value of integrated magnetospheric capacitance to use

!> type of diffusion solver to sue
integer :: diffsolvetype=2       ! 1 - first order backward Euler time stepping; 2 - 2nd order TRBDF2 diffusion solver

!> milestone output information (default to none, i.e. zero value)
integer :: mcadence=-1      ! value less than zero switches this off, > zero gives the cadence at which to perform milestone outputs (in terms of number of outputs per milestone)

!> gravitational drift terms
logical :: flaggravdrift=.false.

!> flag for lagrangian grid (assume drifting at E x B/B**2)
logical :: flaglagrangian=.false.

!> do we consider pressure terms in perp momentum equations
logical :: flagdiamagnetic=.false.

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
!! READS THE INPUT CONFIGURATION FILE, ASSIGNS VARIABLES FOR FILENAMES, SIZES, ETC.

class(gemini_cfg), intent(inout) :: cfg
logical, intent(in), optional :: verbose

integer :: i, realbits, lxp, lyp
real(wp) :: NaN, glat, glon, xdist, ydist, alt_min, alt_max, alt_scale(4), Bincl, nmf, nme

character(:), allocatable :: suffix

!> READ CONFIG FILE FOR THIS SIMULATION
!! NOTE: Namelist file groups must be read in order they appear in the Namelist file, or End of File error occurs
suffix = get_suffix(cfg%infile)

select case (suffix)
case ('.nml')
  call read_nml(cfg, verbose)
case ('.ini')
  call read_ini(cfg)
case default
  error stop 'not sure how to read config file ' // cfg%infile
end select

end subroutine read_configfile


character(5) function get_compiler_vendor() result(vendor)

character(80) :: cvers
integer :: i, j
character(*), parameter :: vendors(2) = [character(5) :: "Intel", "GCC"]

cvers = compiler_version()

do j = 1,size(vendors)
  vendor = vendors(j)
  i = index(cvers, vendor)
  if (i > 0) exit
end do

if(vendor=="GCC") then
  vendor = "GNU"
elseif(i == 0) then
  vendor = ""
  write(stderr,'(A,/,A)') "could not determine compiler vendor from",cvers
end if

end function get_compiler_vendor


function expand_envvar(path, envvar) result(expanded)
!! replace @...@ string like metabuild system e.g. CMake, based on env var
!! this function REQUIRES:
!!
!! 1. envvar is a defined environment variable
!! 2. path contains a matching @envvar@ substring

character(:), allocatable :: expanded, substr
character(*), intent(in) :: path, envvar

integer :: i, L, istat
character(1000) :: buf

expanded = expanduser(path)
!! in case no @var@, and to sanitize fixed width config namelist

i = index(path, "@")
if (i < 1) return
!! nothing to expand

substr = "@" // envvar // "@"
i = index(path, substr)
if (i < 1) return
!! this envvar is not in path, perhaps multiple calls to expand_envvar.
!! Should this be an error?

call get_environment_variable(envvar, buf, length=L, status=istat)
if (istat /= 0 .or. L < 1) error stop "config:expand_envvar: environment variable empty or not defined: " // envvar

expanded = path(1:i - 1) // trim(adjustl(buf)) // path(i + len(substr):len(expanded))

end function expand_envvar

end module config

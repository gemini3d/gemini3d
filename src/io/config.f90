module gemini3d_config

use filesystem, only: suffix
use phys_consts, only : wp

implicit none (type, external)
private
public :: read_configfile, gemini_cfg

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
  integer :: interptype = 0                       ! assumptions to be used when interpolating neutrals:  0-cartesian 2D, 1-axisymmetric 2D, 3 - cartesian 3D
  real(wp) :: sourcemlat = 0,sourcemlon = 0         ! source latitude and longitude
  real(wp) :: dtneu = 0                           ! time step between neutral inputs
  real(wp) :: dxn = 0,drhon = 0,dzn = 0               ! dx,dy (or drho),dz for neutral inputs

  !> preciptiation file inputs
  integer :: flagprecfile                       ! whether or not we have precipitation input from a file
  character(:), allocatable :: precdir          ! location of precipitation input data
  real(wp) :: dtprec = 0                          ! time step between precipitation inputs

  !> electric field file inputs
  integer :: flagE0file                         ! whether or not to have electric field file input
  character(:), allocatable :: E0dir            ! location of electric field input data
  real(wp) :: dtE0 = 0                            ! time step between electric field inputs

  !> solar flux inputs (not implemented yet)
  integer :: flagsolfluxfile=0
  character(:), allocatable :: solfluxdir
  real(wp) :: dtsolflux=0

  !> solar flux inputs (not implemented yet)
  integer :: flagneutralBGfile=0
  character(:), allocatable :: neutralBGdir
  real(wp) :: dtneutralBG=0

  !> GLOW parameters
  integer :: flagglow              ! whether or not to use glow to compute impact ionization
  real(wp) :: dtglow, dtglowout    ! time step between GLOW updates and outputs for GLOW emissions

  !> fang parameters
  integer :: flag_fang = 2008       ! configure Fang ionization model
  integer :: diff_num_flux = 0      ! if flag_fang=0, select input differential number flux type
  real(wp) :: kappa = 1e4_wp        ! for diff_num_flux=1, kappa distribution for kappa > 2
  real(wp) :: bimax_frac = 1._wp    ! for diff_num_flux=2, bimaxwellian with second char. energy bimax_frac * E0
  real(wp) :: W0_char = 3000._wp    ! for diff_num_flux=3, thermal/characteristic energy in eV

  !! parameters below this line can only be changed via the .nml input format
  !> equatorial ionization anomaly
  logical :: flagEIA = .false.       ! whether or not to include and equatorial ionization anomaly in simulation
  real(wp) :: v0equator = 10._wp      ! max vertical drift of plasma at equator for EIA

  !> varying neutral atmosphere background
  logical :: flagneuBG = .false.                ! whether or not to allow MSIS to be called to update neutral background
  real(wp) :: dtneuBG = 900._wp                  ! approximate time between MSIS calls
  integer :: msis_version

  !> background preciptation
  real(wp) :: PhiWBG = 1e-3_wp                      ! background total energy flux in mW/m^2
  real(wp) :: W0BG = 3e3_wp                      ! background characteristic energy for precipitation

  !> parallel current calculations
  logical :: flagJpar = .true.                  ! whether or not to compute parallel current (some simulation setups will give really poor results); code ignores this if potential is resolved along the field line since computing Jpar will not be prone to artifacts as it is in th EFL cases...

  !> inertial capacitance
  integer :: flagcap = 0           ! use inertial capacitance? 0 - set all to zero, 1 - use ionosphere to compute, 2 - add a magnetospheric part
  real(wp) :: magcap = 5        ! value of integrated magnetospheric capacitance to use

  !> type of diffusion solver to sue
  integer :: diffsolvetype = 2       ! 1 - first order backward Euler time stepping; 2 - 2nd order TRBDF2 diffusion solver

  !> milestone output information (default to none, i.e. zero value)
  integer :: mcadence = -1      ! value less than zero switches this off, > zero gives the cadence at which to perform milestone outputs (in terms of number of outputs per milestone)

  !> gravitational drift terms
  logical :: flaggravdrift = .false.

  !> flag for lagrangian grid (assume drifting at E x B/B**2)
  logical :: flaglagrangian = .false.

  !> do we consider pressure terms in perp momentum equations
  logical :: flagdiamagnetic = .false.

  !> do we compute energy and momentum input into the neutrals
  logical :: flagtwoway = .false.

  !> is the background current assumed to be divergence free?
  logical :: flagnodivJ0 = .false.

  !> Farley-Buneman instability
  integer :: flagFBI = 0
  !! default: 0, which does not run FBI model.
  !! 1: turn on only abnormal heating
  !! 2: abnormal heating and non-linear current

  !> electron rotational and vibrational cooling
  integer :: flagevibcool = 0
  !! 1: use new model
  !! 0: use old model

  !> whether or not to compute magnetic pole location based on year (.false.=use default; .true.=compute based on year)
  logical :: flagmagpole=.false.
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

    !> READ CONFIG FILE FOR THIS SIMULATION
    !! NOTE: Namelist file groups must be read in order they appear in the Namelist file, or End of File error occurs
    select case (suffix(cfg%infile))
    case ('.nml')
      call read_nml(cfg, verbose)
    case ('.ini')
      print '(a)', "WARNING: .ini files are long-deprecated and may not work."
      call read_ini(cfg)
    case default
      error stop 'ERROR:gemini3d:config: not sure how to read Gemini3D configuration file: ' // cfg%infile
    end select
  end subroutine read_configfile

end module gemini3d_config

! Audit modification 2026-09-16: deterministic namelist defaults, robust headers and required-value checks.
! Audit modification 2026-09-16: preserve the positive real64 density floor
submodule(gemini3d_config) config_nml

use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use gemini3d_sysinfo, only : expand_envvar, get_compiler_vendor
use filesystem, only : absolute
use phys_consts, only: mindens, mindensnull, mindensdiv
use timeutils, only: ymd2doy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

implicit none (type, external)

contains
  !FIXME:  some default value redundancies below...
  module procedure read_nml
    !! Reads simulation configuration file in .nml
    !! Note that it is best to rewind the file before any read operation, otherwise if the file pointer is already
    !! past the group of interest it will (may?) miss that group and return junk.

    integer :: u, i
    logical :: allow_missing_spatial

    integer :: ymd(3)
    real(wp) :: UTsec0
    real(wp) :: tdur
    real(wp) :: dtout
    real(wp) :: activ(3)
    real(wp) :: tcfl
    real(wp) :: Teinf
    integer :: potsolve, flagperiodic=0, flagoutput, flagcap=0, flag_fang, flagdneu
    integer :: interptype
    real(wp) :: sourcemlat,sourcemlon
    real(wp) :: dtneu
    real(wp) :: dxn=0.0,drhon=0.0,dzn=0.0
    real(wp) :: dtprec=0
    character(1000) :: indat_size, indat_grid, indat_file, source_dir, prec_dir, E0_dir, solfluxdir, neutralBGdir
    character(4) :: file_format=""  !< need to initialize blank or random invisible fouls len_trim>0
    real(wp) :: dtE0=0
    real(wp) :: dtglow=0, dtglowout=0
    logical :: flagEIA
    real(wp) :: v0equator
    real(wp) :: dtsolflux
    real(wp) :: dtneuBGfile

    ! for "default backgroud"
    logical :: flagneuBG=.false.
    real(wp) :: dtneuBG
    integer :: msis_version

    real(wp) :: PhiWBG,W0BG
    logical :: flagJpar
    real(wp) :: magcap
    integer :: diffsolvetype
    integer :: mcadence
    logical :: flaggravdrift
    logical :: flaglagrangian
    logical :: flagdiamagnetic
    logical :: flagtwoway
    logical :: flagnodivJ0

    ! for controlling energy distribution of incident electron flux
    integer :: diff_num_flux
    real(wp) :: kappa, bimax_frac, W0_char

    ! for controlling inclusion of Farley-Buneman anomalous heating/conductance
    integer :: flagFBI

    ! for controlling which electron cooling rates are used in energy equations
    integer :: flagevibcool

    ! user flag to enable calculation of magnetic pole based on year
    logical :: flagmagpole=.false.

    ! controls type of electron velocity solve
    logical :: flagJ1ve

    ! add nightime ionization
    logical :: flagnightQ = .false.

    ! in case the user wants to specify minimum allowed density
    real(wp) :: mindens_userval=1.0e-100_wp
    real(wp) :: mindensnull_userval=1.0e-20
    real(wp) :: mindensdiv_userval=1.0e-5

    namelist /base/ ymd, UTsec0, tdur, dtout, activ, tcfl, Teinf
    namelist /files/ file_format, indat_size, indat_grid, indat_file
    namelist /flags/ potsolve, flagperiodic, flagoutput
    namelist /neutral_perturb/ flagdneu, interptype, sourcemlat, sourcemlon, dtneu, dxn, drhon, dzn, source_dir
    namelist /precip/ dtprec, prec_dir
    namelist /efield/ dtE0, E0_dir
    namelist /fang/ flag_fang
    namelist /glow/ dtglow, dtglowout
    namelist /EIA/ flagEIA,v0equator
    namelist /neutral_BG/ flagneuBG,dtneuBG, msis_version
    namelist /precip_BG/ PhiWBG,W0BG
    namelist /Jpar/ flagJpar
    namelist /capacitance/ flagcap,magcap     ! later need to regroup these in a way that is more logical now there are so many more inputs
    namelist /diffusion/ diffsolvetype
    namelist /milestone/ mcadence
    namelist /gravdrift/ flaggravdrift
    namelist /lagrangian/ flaglagrangian
    namelist /diamagnetic/ flagdiamagnetic
    namelist /twoway_coupled/ flagtwoway
    namelist /nodivJ0/ flagnodivJ0
    namelist /solflux/ dtsolflux,solfluxdir
    namelist /neutralBG_file/ dtneuBGfile, neutralBGdir
    namelist /fang_pars/ diff_num_flux, kappa, bimax_frac, W0_char
    namelist /FBI/ flagFBI
    namelist /evibcool/ flagevibcool
    namelist /magpole/ flagmagpole
    namelist /J1ve/ flagJ1ve
    namelist /nightQ/ flagnightQ
    namelist /input_coverage/ allow_missing_spatial
    namelist /mindens_user/ mindens_userval, mindensnull_userval, mindensdiv_userval

    allow_missing_spatial=.false.

    ! Initialize on every read: declaration initializers have SAVE semantics.
    ymd=0; UTsec0=-1; tdur=-1; dtout=-1; activ=-1; tcfl=-1; Teinf=-1
    potsolve=-1; flagperiodic=0; flagoutput=-1; flagcap=0; flag_fang=2008; flagdneu=1
    interptype=0; sourcemlat=0; sourcemlon=0; dtneu=0; dxn=0; drhon=0; dzn=0
    dtprec=0; dtE0=0; dtglow=0; dtglowout=0; dtsolflux=0; dtneuBGfile=0
    indat_size=""; indat_grid=""; indat_file=""; file_format=""
    source_dir=""; prec_dir=""; E0_dir=""; solfluxdir=""; neutralBGdir=""
    flagEIA=.false.; v0equator=10._wp; flagneuBG=.false.; dtneuBG=900._wp; msis_version=0
    PhiWBG=1e-3_wp; W0BG=3000._wp; flagJpar=.true.; magcap=5._wp
    diffsolvetype=2; mcadence=-1; flaggravdrift=.false.; flaglagrangian=.false.
    flagdiamagnetic=.false.; flagtwoway=.false.; flagnodivJ0=.false.
    diff_num_flux=0; kappa=1e4_wp; bimax_frac=1._wp; W0_char=3000._wp
    flagFBI=0; flagevibcool=0; flagmagpole=.false.; flagJ1ve=.false.; flagnightQ=.false.
    mindens_userval=1e-100_wp; mindensnull_userval=1e-20_wp; mindensdiv_userval=1e-5_wp

    if(.not. allocated(cfg%outdir)) error stop 'gemini3d:config:config_nml please specify simulation output directory'
    if(.not. allocated(cfg%infile)) error stop 'gemini3d:config:config_nml please specify simulation configuration file config.nml'

    open(newunit=u, file=cfg%infile, status='old', action='read')

    read(u, nml=base, iostat=i)
    call check_nml_io(i, cfg%infile, "base")
    if (any(ymd <= 0)) error stop "config: base must specify a positive ymd"
    if (ymd2doy(ymd(1),ymd(2),ymd(3))<1) error stop "config: invalid calendar date"
    if (.not. all(ieee_is_finite([UTsec0,tdur,dtout,activ,tcfl,Teinf]))) &
      error stop "config: nonfinite base value"
    if (UTsec0 < 0 .or. UTsec0 >= 86400 .or. tdur <= 0 .or. dtout <= 0) &
      error stop "config: require 0<=UTsec0<86400 and positive tdur/dtout"
    if (tcfl <= 0 .or. tcfl > 1 .or. Teinf <= 0 .or. any(activ < 0)) &
      error stop "config: invalid/missing tcfl, Teinf or activity indices"
    cfg%ymd0 = ymd
    cfg%UTsec0 = UTsec0
    cfg%tdur = tdur
    cfg%dtout = dtout
    cfg%activ = activ
    cfg%tcfl = tcfl
    cfg%Teinf = Teinf

    rewind(u)
    read(u, nml=flags, iostat=i)
    call check_nml_io(i, cfg%infile, "flags")
    if (.not. any(potsolve == [0,1,3])) &
      error stop "config: potsolve must be 0,1,3; inductive mode 2 is not implemented"
    if (flagoutput < 1 .or. flagoutput > 3) error stop "config: flagoutput must be 1,2,3"
    cfg%potsolve = potsolve
    cfg%flagperiodic = flagperiodic
    cfg%flagoutput = flagoutput

    rewind(u)
    read(u, nml=files, iostat=i)
    call check_nml_io(i, cfg%infile, "files")

    if (len_trim(indat_size)==0 .or. len_trim(indat_grid)==0 .or. len_trim(indat_file)==0) &
      error stop "config: files must specify indat_size, indat_grid, indat_file"

    !> auto file_format if not specified
    if (len_trim(file_format) > 0) then
      cfg%out_format = trim(file_format)
    else
      file_format = suffix(indat_size)
      cfg%out_format = file_format(2:)
    endif

    !> absolute paths or paths relative to cfg%outdir
    ! print '(a)', "TRACE: indat_size " // expand_envvar(indat_size)
    ! print '(a)', "TRACE: outdir = " // cfg%outdir
    cfg%indatsize = absolute(expand_envvar(indat_size), cfg%outdir)
    cfg%indatgrid = absolute(expand_envvar(indat_grid), cfg%outdir)
    cfg%indatfile = absolute(expand_envvar(indat_file), cfg%outdir)
    ! print '(a)', "TRACE: absolute(indat_size) " // cfg%indatsize

    if (namelist_exists(u, "neutral_perturb", verbose)) then
      cfg%flagdneu = 1
      rewind(u)
      read(u, nml=neutral_perturb, iostat=i)
      call check_nml_io(i, cfg%infile, "neutral_perturb")
      if (len_trim(source_dir)==0) error stop "config: missing source_dir"
      cfg%sourcedir = absolute(expand_envvar(source_dir), cfg%outdir)
      cfg%interptype = interptype
      cfg%sourcemlat = sourcemlat
      cfg%sourcemlon = sourcemlon
      cfg%dtneu = dtneu
      cfg%drhon = drhon
      cfg%dzn = dzn
      cfg%dxn = dxn
    else
      cfg%flagdneu = 0
      cfg%sourcedir = ""
    endif

    if (namelist_exists(u, "precip", verbose)) then
      cfg%flagprecfile = 1
      rewind(u)
      read(u, nml=precip, iostat=i)
      call check_nml_io(i, cfg%infile, "precip")
      if (len_trim(prec_dir)==0) error stop "config: missing prec_dir"
      cfg%precdir = absolute(expand_envvar(prec_dir), cfg%outdir)
      cfg%dtprec = dtprec
    else
      cfg%flagprecfile = 0
      cfg%precdir = ""
    endif

    if (namelist_exists(u, "efield", verbose)) then
      cfg%flagE0file = 1
      rewind(u)
      read(u, nml=efield, iostat=i)
      call check_nml_io(i, cfg%infile, "efield")
      if (len_trim(E0_dir)==0) error stop "config: missing E0_dir"
      cfg%E0dir = absolute(expand_envvar(E0_dir), cfg%outdir)
      cfg%dtE0 = dtE0
    else
      cfg%flagE0file = 0
      cfg%E0dir = ""
    endif

    if (namelist_exists(u, "solflux", verbose)) then
      cfg%flagsolfluxfile = 1
      rewind(u)
      read(u, nml=solflux, iostat=i)
      call check_nml_io(i, cfg%infile, "solflux")
      if (len_trim(solfluxdir)==0) error stop "config: missing solfluxdir"
      cfg%solfluxdir = absolute(expand_envvar(solfluxdir), cfg%outdir)
      cfg%dtsolflux = dtsolflux
    else
      cfg%flagsolfluxfile = 0
      cfg%solfluxdir = ""
    endif

    !> neural background (optional)
    if (namelist_exists(u,'neutral_BG')) then
      rewind(u)
      read(u, nml=neutral_BG, iostat=i)
      call check_nml_io(i, cfg%infile, "neutral_BG")
      cfg%flagneuBG=flagneuBG
      cfg%dtneuBG=dtneuBG
      cfg%msis_version = msis_version
    else
      cfg%flagneuBG=.false.
      cfg%msis_version = 0
    end if

    if (namelist_exists(u, "neutralBG_file", verbose)) then
      cfg%flagneutralBGfile = 1
      rewind(u)
      read(u, nml=neutralBG_file, iostat=i)
      call check_nml_io(i, cfg%infile, "neutralBG_file")
      if (len_trim(neutralBGdir)==0) error stop "config: missing neutralBGdir"
      cfg%neutralBGdir = absolute(expand_envvar(neutralBGdir), cfg%outdir)
      cfg%dtneuBGfile = dtneuBGfile
    else
      cfg%flagneutralBGfile = 0
      cfg%neutralBGdir = ""
    endif

    if (namelist_exists(u, "glow", verbose)) then
      cfg%flagglow = 1
      rewind(u)
      read(u, nml=glow, iostat=i)
      call check_nml_io(i, cfg%infile, "glow")
      cfg%dtglow = dtglow
      cfg%dtglowout = dtglowout
    else
      cfg%flagglow = 0
    endif

    !> EIA (optional)
    if (namelist_exists(u,'EIA')) then
      rewind(u)
      read(u, nml=EIA, iostat=i)
      call check_nml_io(i, cfg%infile, "EIA")
      cfg%flagEIA=flagEIA
      cfg%v0equator=v0equator
    else
      cfg%flagEIA=.false.
    end if

    !> precip background (optional)
    if (namelist_exists(u,'precip_BG')) then
      rewind(u)
      read(u, nml=precip_BG, iostat=i)
      call check_nml_io(i, cfg%infile, "precip_BG")
      cfg%PhiWBG=PhiWBG
      cfg%W0BG=W0BG
    else
      cfg%PhiWBG=1e-3_wp
      cfg%W0BG=3000
    end if

    !> parallel current density (optional)
    if (namelist_exists(u,'Jpar')) then
      rewind(u)
      read(u, nml=Jpar, iostat=i)
      call check_nml_io(i, cfg%infile, "Jpar")
      cfg%flagJpar=flagJpar
    else
      cfg%flagJpar=.true.
    end if

    !> inertial capacitance (optional)
    if (namelist_exists(u,'capacitance')) then
      rewind(u)
      read(u, nml=capacitance, iostat=i)
      call check_nml_io(i, cfg%infile, "capacitance")
      cfg%flagcap=flagcap
      cfg%magcap=magcap
    else
      cfg%flagcap=0    !default to zero capacitance
    end if

    !> diffusion solve type (optional). i.e. to switch between backward Euler and TRBDF2
    if (namelist_exists(u,'diffusion')) then
      rewind(u)
      read(u, nml=diffusion, iostat=i)
      call check_nml_io(i, cfg%infile, "diffusion")
      cfg%diffsolvetype=diffsolvetype
    else
      cfg%diffsolvetype=2     !default to TRBDF2 - it almost always works
    end if

    !> information about milestone outputs (optional)
    if (namelist_exists(u,'milestone')) then
      rewind(u)
      read(u,nml=milestone,iostat=i)
      call check_nml_io(i,cfg%infile,"milestone")
      cfg%mcadence = mcadence
    else
      cfg%mcadence = -1     !default to no milestones (<0 is a sentinel value)
    end if

    !> whether or not to include gravitational terms in drift and potential source equations
    if (namelist_exists(u,'gravdrift')) then
      rewind(u)
      read(u,nml=gravdrift,iostat=i)
      call check_nml_io(i,cfg%infile,"gravdrift")
      cfg%flaggravdrift=flaggravdrift
    else
      cfg%flaggravdrift=.false.     !by default do not include grav currents and drifts
    end if

    !> whether or not to allow the grid to drift at the ExB speed
    if (namelist_exists(u,'lagrangian')) then
      rewind(u)
      read(u,nml=lagrangian,iostat=i)
      call check_nml_io(i,cfg%infile,"lagrangian")
      cfg%flaglagrangian=flaglagrangian
    else
      cfg%flaglagrangian=.false.
    end if

    !> whether or not to use pressure terms in perp momentum
    if (namelist_exists(u,'diamagnetic')) then
      rewind(u)
      read(u,nml=diamagnetic,iostat=i)
      call check_nml_io(i,cfg%infile,"diamagnetic")
      cfg%flagdiamagnetic=flagdiamagnetic
    else
      cfg%flagdiamagnetic=.false.
    end if

    !> two-way coupled option
    if (namelist_exists(u,'twoway_coupled')) then
      rewind(u)
      read(u,nml=twoway_coupled,iostat=i)
      call check_nml_io(i,cfg%infile,"twoway_coupled")
      cfg%flagtwoway=flagtwoway
    else
      cfg%flagtwoway=.false.
    end if

    if (namelist_exists(u,'nodivJ0')) then
      rewind(u)
      read(u,nml=nodivJ0,iostat=i)
      call check_nml_io(i,cfg%infile,"nodivJ0")
      cfg%flagnodivJ0=flagnodivJ0
    else
      cfg%flagnodivJ0=.false.
    end if

    if (namelist_exists(u, 'FBI')) then
      rewind(u)
      read(u, nml=FBI, iostat=i)
      call check_nml_io(i, cfg%infile, "FBI")
      cfg%flagFBI = flagFBI
    else
      cfg%flagFBI = 0
    endif

    if (namelist_exists(u, 'evibcool')) then
      rewind(u)
      read(u, nml=evibcool, iostat=i)
      call check_nml_io(i, cfg%infile, "evibcool")
      cfg%flagevibcool = flagevibcool
    else
      cfg%flagevibcool = 0    ! default to legacy rates, for now, so CI still works okay
    endif

    if (namelist_exists(u, 'magpole')) then
      rewind(u)
      read(u, nml=magpole, iostat=i)
      call check_nml_io(i, cfg%infile, "magpole")
      cfg%flagmagpole = flagmagpole
    else
      cfg%flagmagpole = .false.    ! by default use the legacy GEMINI value
    end if

    if (namelist_exists(u, 'J1ve')) then
      rewind(u)
      read(u, nml=J1ve, iostat=i)
      call check_nml_io(i, cfg%infile, "J1ve")
      cfg%flagJ1ve = flagJ1ve
    else
      cfg%flagJ1ve = .false.    ! not incorporating current density into electron drift so CI still works okay
    endif

    if (namelist_exists(u, 'nightQ')) then
      rewind(u)
      read(u, nml=nightQ, iostat=i)
      call check_nml_io(i, cfg%infile, "nightQ")
      cfg%flagnightQ = flagnightQ
    else
      cfg%flagnightQ = .false.    ! not adding nighttime ionization (default uses the older version)
    end if

    if (namelist_exists(u, "fang", verbose)) then
      rewind(u)
      read(u, nml=fang, iostat=i)
      call check_nml_io(i, cfg%infile, "fang")
      cfg%flag_fang = flag_fang
    else
      cfg%flag_fang = 2008  !< legacy default
    endif

    if (namelist_exists(u, "fang_pars", verbose)) then
      rewind(u)
      read(u, nml=fang_pars, iostat=i)
      call check_nml_io(i, cfg%infile, "fang_pars")
      cfg%flag_fang = 0 ! force fang flag for integrated spectrum
      cfg%diff_num_flux = diff_num_flux
      cfg%kappa = kappa
      cfg%bimax_frac = bimax_frac
      cfg%W0_char = W0_char
    else
      cfg%diff_num_flux = 0 ! Maxwellian, same as Fang et al. 2008 within 5% in most cases
      cfg%kappa = 1e4_wp ! close to Maxwellian
      cfg%bimax_frac = 1._wp ! Maxwellian
      cfg%W0_char = 3000._wp ! same as W0BG default
    endif

    if (namelist_exists(u, 'input_coverage')) then
      rewind(u)
      read(u,nml=input_coverage,iostat=i)
      call check_nml_io(i,cfg%infile,'input_coverage')
    endif
    cfg%allow_missing_spatial=allow_missing_spatial

    if (namelist_exists(u, 'mindens_user')) then
      rewind(u)
      read(u, nml=mindens_user, iostat=i)
      call check_nml_io(i, cfg%infile, "mindens_user")
      mindens = mindens_userval    ! this is different from the others since we just directly set the module variable, rather than cfg
      mindensnull = mindensnull_userval
      mindensdiv = mindensdiv_userval
    else
      mindens = 1.0e-100_wp
      mindensnull = 1.0e-20_wp
      mindensdiv  = 1.0e-5_wp
    end if

    if (cfg%flagperiodic<0 .or. cfg%flagperiodic>1) error stop "config: flagperiodic must be 0 or 1"
    if (.not.any(cfg%msis_version==[0,21])) error stop "config: msis_version must be 0 or 21"
    if (.not.any(cfg%flag_fang==[0,2008,2010])) error stop "config: unsupported Fang model"
    if (.not.any(cfg%diffsolvetype==[1,2])) error stop "config: diffusion type must be 1 or 2"
    if (cfg%flagcap<0 .or. cfg%flagcap>2) error stop "config: flagcap must be 0,1,2"
    if (cfg%flagFBI<0 .or. cfg%flagFBI>2) error stop "config: flagFBI must be 0,1,2"
    if (.not.any(cfg%flagevibcool==[0,1])) error stop "config: flagevibcool must be 0 or 1"
    if (cfg%mcadence==0) error stop "config: mcadence must be negative (disabled) or positive"
    if (.not.all(ieee_is_finite([cfg%magcap,cfg%v0equator,cfg%PhiWBG,cfg%W0BG, &
        cfg%kappa,cfg%bimax_frac,cfg%W0_char]))) error stop "config: nonfinite optional value"
    if (cfg%magcap<0 .or. cfg%PhiWBG<0 .or. cfg%W0BG<=0) error stop "config: invalid background/capacitance"
    if (cfg%diff_num_flux<0 .or. cfg%diff_num_flux>4) error stop "config: unsupported spectral distribution"
    if (cfg%kappa<=2 .or. cfg%bimax_frac<=0 .or. cfg%W0_char<=0) error stop "config: invalid Fang parameters"
    if (cfg%flagdneu/=0) then
      if (.not.any(cfg%interptype==[0,1,3,4,5,6])) error stop "config: unsupported neutral interpolation"
      if (.not.all(ieee_is_finite([cfg%sourcemlat,cfg%sourcemlon,cfg%dxn,cfg%drhon,cfg%dzn]))) &
        error stop "config: nonfinite neutral geometry"
      if (abs(cfg%sourcemlat)>90) error stop "config: neutral latitude outside [-90,90]"
    endif
    if (cfg%flagprecfile/=0) call positive_cadence(cfg%dtprec, "dtprec")
    if (cfg%flagE0file/=0) call positive_cadence(cfg%dtE0, "dtE0")
    if (cfg%flagdneu/=0) call positive_cadence(cfg%dtneu, "dtneu")
    if (cfg%flagsolfluxfile/=0) call positive_cadence(cfg%dtsolflux, "dtsolflux")
    if (cfg%flagneutralBGfile/=0) call positive_cadence(cfg%dtneuBGfile, "dtneuBGfile")
    if (cfg%flagneuBG) call positive_cadence(cfg%dtneuBG, "dtneuBG")
    if (cfg%flagglow/=0) then
      call positive_cadence(cfg%dtglow, "dtglow")
      call positive_cadence(cfg%dtglowout, "dtglowout")
    endif
    if (.not.all(ieee_is_finite([mindens,mindensnull,mindensdiv]))) &
      error stop "config: density floors must be finite"
    if (min(mindens,mindensnull,mindensdiv)<=0) error stop "config: density floors must be positive"
    close(u)
  end procedure read_nml

  subroutine positive_cadence(value, name)
    real(wp), intent(in) :: value
    character(*), intent(in) :: name
    if (.not.ieee_is_finite(value)) error stop "config: nonfinite cadence: " // name
    if (value<=0 .or. value>86400) error stop "config: cadence must be in (0,86400]: " // name
  end subroutine positive_cadence


  logical function namelist_exists(u, nml, verbose)
    !! determines if Namelist exists in file

    character(*), intent(in) :: nml    ! FIXME:  is it bad to use a keyword as a variable name?
    integer, intent(in) :: u
    logical, intent(in), optional :: verbose

    logical :: debug
    integer :: i, j, code, n
    character(:), allocatable :: token
    character(256) :: line  !< arbitrary length

    debug = .false.
    if(present(verbose)) debug = verbose

    namelist_exists = .false.

    rewind(u)

    do
      read(u, '(A)', iostat=i) line
      if(i/=0) exit
      line = adjustl(line)
      if (line(1:1) /= '&') cycle
      n = scan(line(2:), " " // achar(9) // "/,!=" )
      if (n==0) n=len_trim(line)
      token = line(2:n)
      do j=1,len(token)
        code=iachar(token(j:j))
        if (code>=iachar('A') .and. code<=iachar('Z')) token(j:j)=achar(code+32)
      enddo
      block
        character(len(nml)) :: expected
        expected=nml
        do j=1,len(expected)
          code=iachar(expected(j:j))
          if (code>=iachar('A') .and. code<=iachar('Z')) expected(j:j)=achar(code+32)
        enddo
        if (token == expected) then
          namelist_exists = .true.
          exit
        endif
      end block
    end do
    rewind(u)

    if (debug) print *, 'namelist ', nml, namelist_exists
  end function namelist_exists


  subroutine check_nml_io(i, filename, namelist)
    !! checks for EOF and gives helpful error
    !! this accommodates non-Fortran 2018 error stop with variable character

    integer, intent(in) :: i
    character(*), intent(in) :: filename
    character(*), intent(in), optional :: namelist
    character(:), allocatable :: nml, msg

    if(i==0) return

    nml = ""
    if(present(namelist)) nml = namelist

    if (is_iostat_end(i)) error stop "namelist " // nml // ': ensure there is a trailing blank line in ' // filename

    msg = ""
    select case (get_compiler_vendor())
    case ("Intel")
      !! https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-list-of-run-time-error-messages
      select case (i)
      case (19)
        msg = "mismatch between variable names in namelist and Fortran code, or problem in variable specification in file"
      case (623)
        msg = "variable specified in Fortran code missing from Namelist file"
      case (17,18,624,625,626,627,628,680,750,759)
        msg = "namelist file format problem"
      end select
    case ("GCC", "GNU")
      select case (i)
      case (5010)
        msg = "mismatch between variable names in namelist and Fortran code, or problem in variable specification in file"
      end select
    end select


    if (len(msg)==0) write(stderr,*) "namelist read error code",i

    error stop 'namelist ' // nml // " from " // filename // " problem: " // msg
  end subroutine check_nml_io
end submodule config_nml

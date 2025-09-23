submodule(gemini3d_config) config_nml

use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use gemini3d_sysinfo, only : expand_envvar, get_compiler_vendor
use filesystem, only : absolute

implicit none (type, external)

contains
  !FIXME:  some default value redundancies below...
  module procedure read_nml
    !! Reads simulation configuration file in .nml
    !! Note that it is best to rewind the file before any read operation, otherwise if the file pointer is already
    !! past the group of interest it will (may?) miss that group and return junk.

    integer :: u, i

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
    real(wp) :: dtneutralBGfile

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
    namelist /neutralBG_file/ dtneutralBGfile, neutralBGdir
    namelist /fang_pars/ diff_num_flux, kappa, bimax_frac, W0_char
    namelist /FBI/ flagFBI
    namelist /evibcool/ flagevibcool
    namelist /magpole/ flagmagpole

    if(.not. allocated(cfg%outdir)) error stop 'gemini3d:config:config_nml please specify simulation output directory'
    if(.not. allocated(cfg%infile)) error stop 'gemini3d:config:config_nml please specify simulation configuration file config.nml'

    open(newunit=u, file=cfg%infile, status='old', action='read')

    read(u, nml=base, iostat=i)
    call check_nml_io(i, cfg%infile, "base")
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
    cfg%potsolve = potsolve
    cfg%flagperiodic = flagperiodic
    cfg%flagoutput = flagoutput

    rewind(u)
    read(u, nml=files, iostat=i)
    call check_nml_io(i, cfg%infile, "files")

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
      cfg%solfluxdir = absolute(expand_envvar(solfluxdir), cfg%outdir)
      cfg%dtsolflux = dtsolflux
    else
      cfg%flagsolfluxfile = 0
      cfg%solfluxdir = ""
    endif

    if (namelist_exists(u, "neutralBG_file", verbose)) then
      cfg%flagneutralBGfile = 1
      rewind(u)
      read(u, nml=solflux, iostat=i)
      call check_nml_io(i, cfg%infile, "neutralBG_file")
      cfg%neutralBGdir = absolute(expand_envvar(neutralBGdir), cfg%outdir)
      cfg%dtneutralBG = dtneutralBGfile
    else
      cfg%flagneutralBGfile = 0
      cfg%neutralBGdir = ""
    endif

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
      read(u, nml=evibcool, iostat=i)
      call check_nml_io(i, cfg%infile, "magpole")
      cfg%flagmagpole = flagmagpole
    else
      cfg%flagmagpole = .false.    ! by default use the legacy GEMINI value
    endif

    close(u)
  end procedure read_nml


  logical function namelist_exists(u, nml, verbose)
    !! determines if Namelist exists in file

    character(*), intent(in) :: nml    ! FIXME:  is it bad to use a keyword as a variable name?
    integer, intent(in) :: u
    logical, intent(in), optional :: verbose

    logical :: debug
    integer :: i
    character(256) :: line  !< arbitrary length

    debug = .false.
    if(present(verbose)) debug = verbose

    namelist_exists = .false.

    rewind(u)

    do
      read(u, '(A)', iostat=i) line
      if(i/=0) exit
      if (line(1:1) /= '&') cycle
      if (line(2:) == nml) then
        namelist_exists = .true.
        exit
      end if
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

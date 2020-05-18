module config

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, compiler_version

use pathlib, only : expanduser, get_suffix
use phys_consts, only : wp

implicit none (type, external)
private
public :: read_configfile, gemini_cfg

type :: gemini_cfg

!> basic simulation information
integer, dimension(3) :: ymd0
real(wp) :: UTsec0, tdur, dtout
real(wp) :: activ(3)
real(wp) :: tcfl
real(wp) :: Teinf
integer :: potsolve,flagperiodic,flagoutput
logical :: nooutput = .false.
logical :: dryrun = .false.
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

!> equatorial ionization anomaly
logical :: flagEIA=.false.       ! whether or not to include and equatorial ionization anomaly in simulation
real(wp) :: v0equator=10._wp      ! max vertical drift of plasma at equator for EIA

!> varying neutral atmosphere background
logical :: flagneuBG=.false.                ! whether or not to allow MSIS to be called to update neutral background
real(wp) :: dtneuBG=900._wp                  ! approximate time between MSIS calls

!> background preciptation
real(wp) :: PhiWBG=0.05                      ! background total energy flux in mW/m^2  
real(wp) :: W0BG=3e3_wp                      ! background characteristic energy for precipitation

!> parallel current calculations
logical :: flagJpar=.true.                  ! whether or not to compute parallel current (some simulation setups will give really poor results); code ignores this if potential is resolved along the field line since computing Jpar will not be prone to artifacts as it is in th EFL cases...

!> inertial capacitance
integer :: flagcap               ! use inertial capacitance? 0 - set all to zero, 1 - use ionosphere to compute, 2 - add a magnetospheric part
real(wp) :: magcap=30._wp        ! value of integrated magnetospheric capacitance to use 

end type gemini_cfg

character(:), allocatable :: compiler_vendor

contains

subroutine read_configfile(cfg)

class(gemini_cfg), intent(inout) :: cfg

!! READS THE INPUT CONFIGURAITON FILE, ASSIGNS VARIABLES FOR FILENAMES, SIZES, ETC.
integer :: i, realbits, lxp, lyp
real(wp) :: NaN, glat, glon, xdist, ydist, alt_min, alt_max, alt_scale(4), Bincl, nmf, nme

compiler_vendor = get_compiler_vendor()

!> READ CONFIG FILE FOR THIS SIMULATION
!! NOTE: Namelist file groups must be read in order they appear in the Namelist file, or End of File error occurs

if (cfg%infile(len(cfg%infile)-3 : len(cfg%infile)) == '.nml') then
  call read_nml(cfg)
else
  call read_ini(cfg)
endif

end subroutine read_configfile


subroutine read_nml(cfg)

class(gemini_cfg), intent(inout) :: cfg

integer :: u, i

integer :: ymd(3)
real(wp) :: UTsec0
real(wp) :: tdur
real(wp) :: dtout
real(wp) :: activ(3)
real(wp) :: tcfl
real(wp) :: Teinf
integer :: potsolve, flagperiodic=0, flagoutput, flagcap=0
integer :: flagdneu=0
integer :: interptype
real(wp) :: sourcemlat,sourcemlon
real(wp) :: dtneu
real(wp) :: dxn,drhon,dzn
integer :: flagprecfile=0
real(wp) :: dtprec=0
character(256) :: indat_size, indat_grid, indat_file, source_dir, prec_dir, E0_dir
character(4) :: file_format=""  !< need to initialize blank or random invisible fouls len_trim>0
integer :: flagE0file=0
real(wp) :: dtE0=0
integer :: flagglow=0
real(wp) :: dtglow=0, dtglowout=0

namelist /base/ ymd, UTsec0, tdur, dtout, activ, tcfl, Teinf
namelist /files/ file_format, indat_size, indat_grid, indat_file
namelist /flags/ potsolve, flagperiodic, flagoutput, flagcap, flagdneu, flagprecfile, flagE0file, flagglow
namelist /neutral_perturb/ interptype, sourcemlat, sourcemlon, dtneu, dxn, drhon, dzn, source_dir
namelist /precip/ dtprec, prec_dir
namelist /efield/ dtE0, E0_dir
namelist /glow/ dtglow, dtglowout

open(newunit=u, file=cfg%infile, status='old', action='read')

read(u, nml=base, iostat=i)
call check_nml_io(i, cfg%infile, "base", compiler_vendor)
cfg%ymd0 = ymd
cfg%UTsec0 = UTsec0
cfg%tdur = tdur
cfg%dtout = dtout
cfg%activ = activ
cfg%tcfl = tcfl
cfg%Teinf = Teinf

read(u, nml=flags, iostat=i)
call check_nml_io(i, cfg%infile, "flags", compiler_vendor)
cfg%potsolve = potsolve
cfg%flagperiodic = flagperiodic
cfg%flagoutput = flagoutput
cfg%flagcap = flagcap
cfg%flagdneu = flagdneu
cfg%flagprecfile = flagprecfile
cfg%flagE0file = flagE0file
cfg%flagglow = flagglow

read(u, nml=files, iostat=i)
call check_nml_io(i, cfg%infile, "files", compiler_vendor)

!> auto file_format if not specified
if (len_trim(file_format) > 0) then
  cfg%out_format = trim(file_format)
else
  file_format = get_suffix(indat_size)
  cfg%out_format = file_format(2:)
endif

cfg%indatsize = expanduser(indat_size)
cfg%indatgrid = expanduser(indat_grid)
cfg%indatfile = expanduser(indat_file)

if (cfg%flagdneu == 1) then
  read(u, nml=neutral_perturb, iostat=i)
  call check_nml_io(i, cfg%infile, "neutral_perturb", compiler_vendor)
  cfg%sourcedir = expanduser(source_dir)
  cfg%interptype = interptype
  cfg%sourcemlat = sourcemlat
  cfg%sourcemlon = sourcemlon
  cfg%dtneu = dtneu
  cfg%drhon = drhon
  cfg%dzn = dzn
  cfg%dxn = dxn
else
  cfg%sourcedir = ""
endif


if (cfg%flagprecfile == 1) then
  read(u, nml=precip, iostat=i)
  call check_nml_io(i, cfg%infile, "precip", compiler_vendor)
  cfg%precdir = expanduser(prec_dir)
  cfg%dtprec = dtprec
else
  cfg%precdir = ""
endif

if (cfg%flagE0file == 1) then
  read(u, nml=efield, iostat=i)
  call check_nml_io(i, cfg%infile, "efield", compiler_vendor)
  cfg%E0dir = expanduser(E0_dir)
  cfg%dtE0 = dtE0
else
  cfg%E0dir = ""
endif

if (cfg%flagglow == 1) then
  read(u, nml=glow, iostat=i)
  call check_nml_io(i, cfg%infile, "glow", compiler_vendor)
  cfg%dtglow = dtglow
  cfg%dtglowout = dtglowout
endif


close(u)

end subroutine read_nml


subroutine read_ini(cfg)

class(gemini_cfg), intent(inout) :: cfg

integer :: u,i
character(256) :: buf

open(newunit=u, file=cfg%infile, status='old', action='read')

read(u,*) cfg%ymd0(3), cfg%ymd0(2), cfg%ymd0(1)
read(u,*) cfg%UTsec0
read(u,*) cfg%tdur
read(u,*) cfg%dtout
read(u,*) cfg%activ(1), cfg%activ(2), cfg%activ(3)
read(u,*) cfg%tcfl
read(u,*) cfg%Teinf
read(u,*) cfg%potsolve
read(u,*) cfg%flagperiodic
read(u,*) cfg%flagoutput
read(u,*) cfg%flagcap  !< line 11 config.ini
read(u,'(a256)') buf
!! format specifier needed, else it reads just one character
cfg%indatsize = expanduser(buf)
read(u,'(a256)') buf
cfg%indatgrid = expanduser(buf)
read(u,'(a256)') buf
cfg%indatfile = expanduser(buf)   !< line 14
cfg%out_format = "raw"

!! neutral
read(u,*, iostat=i) cfg%flagdneu  !< line 15
call check_ini_io(i, cfg%infile)
if(cfg%flagdneu==1) then
  read(u,*) cfg%interptype
  read(u,*) cfg%sourcemlat, cfg%sourcemlon
  read(u,*) cfg%dtneu
  if (cfg%interptype==3) then     !< read in extra dxn if 3D
    read(u,*) cfg%dxn,cfg%drhon,cfg%dzn
  else
    read(u,*) cfg%drhon,cfg%dzn
  end if
  read(u,'(A256)') buf
  cfg%sourcedir = expanduser(buf)
else
  cfg%sourcedir = ""
endif

read(u,*, iostat=i) cfg%flagprecfile
call check_ini_io(i, cfg%infile)
if (cfg%flagprecfile==1) then
!! get the location of the precipitation input files
  read(u,*, iostat=i) cfg%dtprec
  read(u,'(A256)', iostat=i) buf
  call check_nml_io(i, cfg%infile)
  cfg%precdir = expanduser(buf)
else
  cfg%precdir = ""
end if

read(u,*, iostat=i) cfg%flagE0file
call check_ini_io(i, cfg%infile)
if (cfg%flagE0file==1) then
!! get the location of the precipitation input files
  read(u,*, iostat=i) cfg%dtE0
  read(u,'(a256)', iostat=i) buf
  call check_ini_io(i, cfg%infile)
  cfg%E0dir = expanduser(buf)
else
  cfg%E0dir = ""
end if

read(u,*, iostat=i) cfg%flagglow
call check_ini_io(i, cfg%infile)
if (cfg%flagglow==1) then
  read(u,*, iostat=i) cfg%dtglow
  read(u,*, iostat=i) cfg%dtglowout
  call check_nml_io(i, cfg%infile)
end if

close(u)

end subroutine read_ini


subroutine check_nml_io(i, filename, group, vendor)
!! checks for EOF and gives helpful error
!! this accomodates non-Fortran 2018 error stop with variable character

integer, intent(in) :: i
character(*), intent(in) :: filename
character(*), intent(in), optional :: group, vendor

character(:), allocatable :: grp, msg

if(i==0) return

grp = ""
if(present(group)) grp = group

if (is_iostat_end(i)) then
  write(stderr,*) 'ERROR: group ' // grp // ': ensure there is a trailing blank line in ' // filename
  write(stderr,*) 'also, ensure that Namelist groups are read in order: base, flags, files etc'
  error stop 5
endif

msg = ""
if (present(vendor)) then
  select case (vendor)
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
  case ("GCC")
    select case (i)
    case (5010)
      msg = "mismatch between variable names in namelist and Fortran code, or problem in variable specification in file"
    end select
  end select
endif

if (len(msg)==0) write(stderr,*) "namelist read error code",i

write(stderr,'(A,/,A)') 'ERROR: reading ' // grp // " " // filename, msg
error stop 5

end subroutine check_nml_io


subroutine check_ini_io(i, filename)
!! checks for EOF and gives helpful error
!! this accomodates non-Fortran 2018 error stop with variable character

integer, intent(in) :: i
character(*), intent(in) :: filename

if (is_iostat_end(i)) return

if (i /= 0) then
write(stderr,*) 'ERROR: problem reading ' // filename
error stop 5
endif

end subroutine check_ini_io


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

if(allocated(vendor)) return

vendor = ""
write(stderr,'(A,/,A)') "could not determine compiler vendor from",cvers

end function get_compiler_vendor

end module config

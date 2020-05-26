module config

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, compiler_version

use pathlib, only : expanduser, get_suffix
use phys_consts, only : wp

implicit none (type, external)
private
public :: read_configfile, gemini_cfg

type :: gemini_cfg

integer :: ymd(3)
real(wp) :: UTsec0, tdur, dtout
real(wp) :: activ(3)
real(wp) :: tcfl
real(wp) :: Teinf
integer :: potsolve, flagperiodic, flagoutput, flagcap,flagdneu,flagprecfile,interptype=0
real(wp) :: sourcemlat=0,sourcemlon=0,dtneu=0, dxn=0,drhon=0,dzn=0
logical :: nooutput = .false.
logical :: dryrun = .false.
real(wp) :: dtprec=0
character(:), allocatable :: infile, outdir, &
  indatsize,indatgrid, indatfile, sourcedir, precdir, E0dir, out_format
integer :: flagE0file
real(wp) :: dtE0=0
integer :: flagglow
real(wp) :: dtglow, dtglowout

end type gemini_cfg


interface
module subroutine read_nml(cfg)
class(gemini_cfg), intent(inout) :: cfg
end subroutine read_nml
end interface


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



subroutine read_ini(cfg)

class(gemini_cfg), intent(inout) :: cfg

integer :: u,i
character(256) :: buf

open(newunit=u, file=cfg%infile, status='old', action='read')

read(u,*) cfg%ymd(3), cfg%ymd(2), cfg%ymd(1)
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
  call check_ini_io(i, cfg%infile)
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
  call check_ini_io(i, cfg%infile)
end if

close(u)

end subroutine read_ini


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

module config

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, compiler_version

use pathlib, only : expanduser, get_suffix
use phys_consts, only : wp

implicit none (type, external)
private
public :: read_configfile, gemini_cfg, get_compiler_vendor

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

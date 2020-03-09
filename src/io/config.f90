module config

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, compiler_version
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan

use pathlib, only : expanduser
use phys_consts, only : wp

implicit none

contains

subroutine read_configfile(infile,ymd,UTsec0,tdur,dtout,activ,tcfl,Teinf, &
  potsolve,flagperiodic, flagoutput,flagcap,&
  indatsize,indatgrid,indatfile,flagdneu,interptype, &
  sourcemlat,sourcemlon,dtneu,dxn,drhon,dzn,sourcedir,flagprecfile,&
  dtprec,precdir,flagE0file,dtE0,E0dir,flagglow,dtglow,dtglowout, out_format)

character(*), intent(in) :: infile
integer, dimension(3), intent(out):: ymd
real(wp), intent(out) :: UTsec0
real(wp), intent(out) :: tdur
real(wp), intent(out) :: dtout
real(wp), dimension(3), intent(out) :: activ
real(wp), intent(out) :: tcfl
real(wp), intent(out) :: Teinf
integer, intent(out) :: potsolve, flagperiodic, flagoutput, flagcap
integer, intent(out) :: flagdneu
integer, intent(out) :: interptype
real(wp), intent(out) :: sourcemlat,sourcemlon
real(wp), intent(out) :: dtneu
real(wp), intent(out) :: dxn,drhon,dzn
integer, intent(out) :: flagprecfile
real(wp), intent(out) :: dtprec
character(:), allocatable, intent(out) :: indatsize,indatgrid, indatfile, sourcedir, precdir, E0dir, out_format
integer, intent(out) :: flagE0file
real(wp), intent(out) :: dtE0
integer, intent(out) :: flagglow
real(wp), intent(out) :: dtglow, dtglowout
!! READS THE INPUT CONFIGURAITON FILE, ASSIGNS VARIABLES FOR FILENAMES, SIZES, ETC.
character(256) :: buf, indat_size, indat_grid, indat_file, source_dir, prec_dir, E0_dir
character(4) :: file_format
integer :: i, realbits, lxp, lyp
real(wp) :: NaN, glat, glon, xdist, ydist, alt_min, alt_max, alt_scale(4), Bincl, nmf, nme
logical :: is_nml
character(:), allocatable :: compiler_vendor

namelist /base/ ymd, UTsec0, tdur, dtout, activ, tcfl, Teinf
namelist /files/ file_format, indat_size, indat_grid, indat_file
namelist /flags/ potsolve, flagperiodic, flagoutput, flagcap, flagdneu, flagprecfile, flagE0file, flagglow
namelist /setup/ glat, glon, xdist, ydist, alt_min, alt_max, alt_scale,lxp,lyp,Bincl,nmf,nme
namelist /neutral_perturb/ interptype, sourcemlat, sourcemlon, dxn, drhon, dzn, source_dir
namelist /precip/ dtprec, prec_dir
namelist /efield/ dtE0, E0_dir
namelist /glow/ dtglow, dtglowout

compiler_vendor = get_compiler_vendor()

NaN = ieee_value(0._wp, ieee_quiet_nan)

is_nml = infile(len(infile)-3:len(infile)) == '.nml'

!> READ CONFIG FILE FOR THIS SIMULATION
rawconfig : block
integer :: u
open(newunit=u, file=infile, status='old', action='read')
if (is_nml) then

  read(u, nml=base, iostat=i)
  call check_nml_io(i, infile, "base", compiler_vendor)

  read(u, nml=files, iostat=i)
  call check_nml_io(i, infile, "files", compiler_vendor)
  out_format = trim(file_format)
  indatsize = expanduser(indat_size)
  indatgrid = expanduser(indat_grid)
  indatfile = expanduser(indat_file)

  read(u, nml=flags, iostat=i)
  call check_nml_io(i, infile, "flags", compiler_vendor)

else
  read(u,*) ymd(3),ymd(2),ymd(1)
  read(u,*) UTsec0
  read(u,*) tdur
  read(u,*) dtout
  read(u,*) activ(1),activ(2),activ(3)
  read(u,*) tcfl
  read(u,*) Teinf
  read(u,*) potsolve
  read(u,*) flagperiodic
  read(u,*) flagoutput
  read(u,*) flagcap  ! line 11 config.ini
  read(u,'(a256)') buf  !! format specifier needed, else it reads just one character
  indatsize = expanduser(buf)
  read(u,'(a256)') buf
  indatgrid = expanduser(buf)
  read(u,'(a256)') buf
  indatfile = expanduser(buf)   ! line 14
  out_format = "raw"
endif

!> NEUTRAL PERTURBATION INPUT INFORMATION
!> defaults
interptype=0
sourcemlat=0
sourcemlon=0
dtneu=0
drhon=0
dzn=0
dxn=0
if(is_nml) then
  if (flagdneu == 1) then
    read(u, nml=neutral_perturb, iostat=i)
    call check_nml_io(i, infile, "neutral_perturb", compiler_vendor)
    sourcedir = expanduser(source_dir)
  endif
else
  read(u,*, iostat=i) flagdneu  ! line 15
  call check_ini_io(i, infile)
  if( flagdneu==1) then
    read(u,*) interptype
    read(u,*) sourcemlat,sourcemlon
    read(u,*) dtneu
    if (interptype==3) then     !read in extra dxn if 3D
      read(u,*) dxn,drhon,dzn
    else
      read(u,*) drhon,dzn
    end if
    read(u,'(A256)') buf
    sourcedir = expanduser(buf)
  endif
end if
!> have to allocate, even when not used, to avoid runtime errors with pickier compilers
if (flagdneu/=1) sourcedir = ""

!> PRECIPITATION FILE INPUT INFORMATION
!> defaults
dtprec=0
if(is_nml) then
  if (flagprecfile == 1) then
    read(u, nml=precip, iostat=i)
    call check_nml_io(i, infile, "precip", compiler_vendor)
    precdir = expanduser(prec_dir)
  endif
else
  read(u,*, iostat=i) flagprecfile
  call check_ini_io(i, infile)
  if (flagprecfile==1) then
  !! get the location of the precipitation input files
    read(u,*, iostat=i) dtprec
    read(u,'(A256)', iostat=i) buf
    call check_nml_io(i, infile)
    precdir = expanduser(buf)
  end if
end if
!> have to allocate, even when not used, to avoid runtime errors with pickier compilers
if (flagprecfile/=1) precdir = ""

!> ELECTRIC FIELD FILE INPUT INFORMATION
!> defaults
dtE0=0
if(is_nml) then
  if (flagE0file == 1) then
    read(u, nml=efield, iostat=i)
    call check_nml_io(i, infile, "efield", compiler_vendor)
    E0dir = expanduser(E0_dir)
  endif
else
  read(u,*, iostat=i) flagE0file
  call check_ini_io(i, infile)
  if (flagE0file==1) then
  !! get the location of the precipitation input files
    read(u,*, iostat=i) dtE0
    read(u,'(a256)', iostat=i) buf
    call check_ini_io(i, infile)
    E0dir = expanduser(buf)
  end if
end if
!> have to allocate, even when not used, to avoid runtime errors with pickier compilers
if (flagE0file/=1) E0dir = ""

!> GLOW ELECTRON TRANSPORT INFORMATION
!> defaults
dtglow=NaN
dtglowout=NaN
if(is_nml) then
  if (flagglow == 1) then
    read(u, nml=glow, iostat=i)
    call check_nml_io(i, infile, "glow", compiler_vendor)
  endif
else
  read(u,*, iostat=i) flagglow
  call check_ini_io(i, infile)
  if (flagglow==1) then
    read(u,*, iostat=i) dtglow
    read(u,*, iostat=i) dtglowout
    call check_nml_io(i, infile)
  end if
end if

close(u)
end block rawconfig

end subroutine read_configfile


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

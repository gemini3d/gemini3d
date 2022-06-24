submodule (gemini3d_config) config_ini

use filesystem, only : expanduser

implicit none (type, external)

contains

module procedure read_ini

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
cfg%out_format = "dat"

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

end procedure read_ini


subroutine check_ini_io(i, filename)
!! checks for EOF and gives helpful error
!! this accommodates non-Fortran 2018 error stop with variable character

integer, intent(in) :: i
character(*), intent(in) :: filename

if (is_iostat_end(i)) return

if (i /= 0) error stop 'ERROR: problem reading ' // filename

end subroutine check_ini_io

end submodule config_ini

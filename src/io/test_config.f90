program test_config
!! test config file reading from Fortran, as NAMELIST has its quirks
!! the order of variables in namelist specification doesn't have to match that of file namelist.

use config, only : read_configfile
use phys_consts, only : wp

implicit none

character(256) :: argv
integer :: i

integer :: ymd(3)
real(wp) :: UTsec0
real(wp) :: tdur
real(wp) :: dtout
real(wp) :: activ(3)
real(wp) :: tcfl
real(wp) :: Teinf
integer :: potsolve, flagperiodic, flagoutput, flagcap
integer :: flagdneu
integer :: interptype
real(wp) :: sourcemlat,sourcemlon
real(wp) :: dtneu
real(wp) :: dxn,drhon,dzn
integer :: flagprecfile
real(wp) :: dtprec
character(:), allocatable :: infile, indatsize,indatgrid, indatfile, sourcedir, precdir, E0dir, out_format
integer :: flagE0file
real(wp) :: dtE0
integer :: flagglow
real(wp) :: dtglow, dtglowout

call get_command_argument(1, argv, status=i)
if (i/=0) error stop 77

infile = trim(argv)

call read_configfile(infile,ymd,UTsec0,tdur,dtout,activ,tcfl,Teinf, &
  potsolve,flagperiodic, flagoutput,flagcap,&
  indatsize,indatgrid,indatfile,flagdneu,interptype, &
  sourcemlat,sourcemlon,dtneu,dxn,drhon,dzn,sourcedir,flagprecfile,&
  dtprec,precdir,flagE0file,dtE0,E0dir,flagglow,dtglow,dtglowout, out_format)

print *, "OK: config read ",infile

end program

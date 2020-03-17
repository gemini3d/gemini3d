program test_config
!! test config file reading from Fortran, as NAMELIST has its quirks
!! the order of variables in namelist specification doesn't have to match that of file namelist.

use config, only : read_configfile, gemini_cfg

implicit none

type(gemini_cfg) :: cfg

character(256) :: argv
character(:), allocatable :: infile
integer :: i

call get_command_argument(1, argv, status=i)
if (i/=0) error stop 77

infile = trim(argv)

call read_configfile(infile,cfg)

print *, "OK: config read ",infile

end program

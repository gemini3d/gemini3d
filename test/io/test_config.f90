program test_config
!! test config file reading from Fortran, as NAMELIST has its quirks
!! the order of variables in namelist specification doesn't have to match that of file namelist.

use gemini3d_config, only : read_configfile, gemini_cfg

implicit none (type, external)

type(gemini_cfg) :: cfg

character(256) :: argv
integer :: i

call get_command_argument(1, argv, status=i)
if (i/=0) error stop "please specify .nml file to read"

cfg%infile = trim(argv)

call read_configfile(cfg, verbose=.true.)

print *, "OK: config read ",cfg%infile

end program

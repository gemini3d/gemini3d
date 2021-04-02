program test_nml
!! test reference config.nml to see if it's a problem with file or Gemini

implicit none (type, external)

call print_nml()
print *, "OK: read basic namelist file"

contains


subroutine print_nml()

type :: gemini_config
  integer :: potsolve, flagperiodic, flagoutput, flagcap
end type gemini_config

character(256) :: argv
integer :: i, u
character(:), allocatable :: infile
type(gemini_config) :: cfg

integer :: potsolve, flagperiodic, flagoutput, flagcap
namelist /flags/ potsolve, flagperiodic, flagoutput, flagcap

call get_command_argument(1, argv, status=i)
if (i/=0) error stop "please specify .nml file to read"
infile = trim(argv)

open(newunit=u, file=infile, status='old', action='read')

read(u, nml=flags)
cfg%potsolve = potsolve
cfg%flagperiodic = flagperiodic
cfg%flagoutput = flagoutput
cfg%flagcap = flagcap
print *,cfg

close(u)

end subroutine print_nml

end program

program test_nml
!! test reference config.nml to see if it's a problem with file or Gemini

implicit none (type, external)

type :: gemini_config
  integer :: potsolve, flagperiodic, flagoutput, flagcap, flagdneu, flagprecfile, flagE0file, flagglow
end type gemini_config

character(256) :: argv
integer :: i, u
character(:), allocatable :: infile
type(gemini_config) :: cfg

integer :: potsolve, flagperiodic, flagoutput, flagcap, flagdneu, flagprecfile, flagE0file, flagglow
namelist /flags/ potsolve, flagperiodic, flagoutput, flagcap, flagdneu, flagprecfile, flagE0file, flagglow

call get_command_argument(1, argv, status=i)
if (i/=0) error stop 77
infile = trim(argv)

open(newunit=u, file=infile, status='old', action='read')

read(u, nml=flags)
cfg%potsolve = potsolve
cfg%flagperiodic = flagperiodic
cfg%flagoutput = flagoutput
cfg%flagcap = flagcap
cfg%flagdneu = flagdneu
cfg%flagprecfile = flagprecfile
cfg%flagE0file = flagE0file
cfg%flagglow = flagglow
print *,cfg

close(u)

end program
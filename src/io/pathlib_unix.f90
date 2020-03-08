submodule (pathlib) pathlib_unix

implicit none

contains

module procedure copyfile

character(len(source)) :: src
character(len(dest)) :: dst
logical :: exists
integer :: icstat
!! https://linux.die.net/man/1/cp
character(*), parameter ::  CMD='cp -rf '
src = source
dst = dest

call execute_command_line(CMD//src//' '//dst, exitstat=istat, cmdstat=icstat)
if (istat == 0 .and. icstat /= 0) istat = icstat
if (istat /= 0) write(stderr,*) 'ERROR: '//CMD//src//' '//dst

end procedure copyfile


module procedure mkdir
!! create a directory, with parents if needed
character(len(path)) :: p
integer :: icstat

character(*), parameter ::  CMD='mkdir -p '
p = path

call execute_command_line(CMD // p, exitstat=istat, cmdstat=icstat)
if (istat == 0 .and. icstat /= 0) istat = icstat
if (istat /= 0) write(stderr,*) 'ERROR: '//CMD//p

end procedure mkdir

end submodule pathlib_unix
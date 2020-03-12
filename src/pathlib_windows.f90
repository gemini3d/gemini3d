submodule (pathlib) pathlib_windows

implicit none

contains

module procedure copyfile

character(len(source)) :: src
character(len(dest)) :: dst
logical :: exists
integer :: icstat
!! https://docs.microsoft.com/en-us/windows-server/administration/windows-commands/copy
character(*), parameter :: CMD='copy /y '

src = filesep_swap(source)
dst = filesep_swap(dest)


call execute_command_line(CMD//src//' '//dst, exitstat=istat, cmdstat=icstat)
if (istat == 0 .and. icstat /= 0) istat = icstat
if (istat /= 0) write(stderr,*) 'ERROR: '//CMD//src//' '//dst

end procedure copyfile


module procedure mkdir
!! create a directory, with parents if needed
character(len(path)) :: p
integer :: icstat
!! https://docs.microsoft.com/en-us/windows-server/administration/windows-commands/md
character(*), parameter :: CMD='mkdir '
p = filesep_swap(path)

call execute_command_line(CMD // p, exitstat=istat, cmdstat=icstat)
if (istat == 0 .and. icstat /= 0) istat = icstat
if (istat /= 0) write(stderr,*) 'ERROR: '//CMD//p

end procedure mkdir


end submodule pathlib_windows
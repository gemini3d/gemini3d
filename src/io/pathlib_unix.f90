submodule (pathlib) pathlib_unix

implicit none (type, external)

contains


module procedure is_absolute

is_absolute = .false.
if(len_trim(path) > 0) is_absolute = path(1:1) == "/"

end procedure is_absolute


module procedure copyfile

integer :: i, j
!! https://linux.die.net/man/1/cp
character(*), parameter ::  CMD='cp -rf '

call execute_command_line(CMD // source // ' ' // dest, exitstat=i, cmdstat=j)
if (i /= 0 .or. j /= 0) error stop "could not copy " // source // " => " // dest

end procedure copyfile


module procedure mkdir
!! create a directory, with parents if needed
integer :: i, j

character(*), parameter ::  CMD='mkdir -p '

if(directory_exists(path)) return

call execute_command_line(CMD // path, exitstat=i, cmdstat=j)
if (i /= 0 .or. j /= 0) error stop "could not create directory " // path

end procedure mkdir

end submodule pathlib_unix

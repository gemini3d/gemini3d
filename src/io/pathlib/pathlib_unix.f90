submodule (pathlib) pathlib_unix

implicit none (type, external)

contains


module procedure is_absolute
!! is path absolute
character(:), allocatable  :: buf

is_absolute = .false.

buf = expanduser(path)

if(len_trim(buf) > 0) is_absolute = buf(1:1) == "/"

end procedure is_absolute


module procedure copyfile
!! copy file from src to dst
!!
!! https://linux.die.net/man/1/cp
integer :: i, j
character(*), parameter ::  CMD='cp -rf '

character(:), allocatable  :: s,d

s = expanduser(source)
d = expanduser(dest)

call execute_command_line(CMD // s // ' ' // d, exitstat=i, cmdstat=j)
if (i /= 0 .or. j /= 0) error stop "could not copy " // source // " => " // dest

end procedure copyfile


module procedure mkdir
!! create a directory, with parents if needed
integer :: i, j

character(*), parameter ::  CMD='mkdir -p '

character(:), allocatable  :: buf

buf = expanduser(path)

if(directory_exists(buf)) return

call execute_command_line(CMD // buf, exitstat=i, cmdstat=j)
if (i /= 0 .or. j /= 0) error stop "could not create directory " // path

end procedure mkdir

end submodule pathlib_unix

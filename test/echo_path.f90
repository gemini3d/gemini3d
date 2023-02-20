program echo_path
!! tell environment variable.
!! can help establish if run frameworks like CTest are passing along variables properly.

implicit none (type, external)

character(:), allocatable :: var
character(4096) :: buf
integer :: i

if (command_argument_count() == 0) then
  var = "PATH"
else
  call get_command_argument(1, buf)
  var = trim(buf)
endif

call get_environment_variable(var, buf, status=i)
if (i==1) error stop var // " is not an environment variable."
if (i==-1) error stop "buffer is too small for " // var
if (i/=0) error stop "problem reading " // var

print *, trim(buf)

end program

program test_expand

use gemini3d_sysinfo, only : expand_envvar

implicit none (type, external)


call test_expand_envvar()
print *, "OK: expand env var"

contains

subroutine test_expand_envvar()

character(*), parameter :: str="abc@test__gem@xyz"
character(11) :: ret

ret = expand_envvar(str)

if (ret /= "abchelloxyz") error stop "expand_envar: expected abchelloxyz, got " // ret

end subroutine test_expand_envvar


end program

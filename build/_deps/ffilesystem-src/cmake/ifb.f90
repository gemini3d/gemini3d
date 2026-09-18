program ifb_check

use, intrinsic :: iso_c_binding

implicit none

interface
integer(C_INT) function realpath_cfi(path, result) bind(C)
import C_INT
character(*), intent(in)  :: path
character(*), intent(out) :: result
end function
end interface

character(:), allocatable :: out
integer(C_INT) :: rc

allocate(character(1000) :: out)

rc = realpath_cfi(".", out)

if (rc < 0) error stop "realpath_cfi failed"

print '(i0)', rc

end program ifb_check

module cpu

use, intrinsic :: iso_c_binding, only : c_int
use, intrinsic :: iso_fortran_env, only : int32
implicit none (type, external)

interface
integer(c_int) function cpu_count_c() bind(c, name="cpu_count")
import c_int
end function cpu_count_c
end interface

contains

integer function cpu_count()

cpu_count = int(cpu_count_c(), int32)

end function cpu_count


end module cpu

submodule (runner) sysdet

use, intrinsic :: iso_c_binding, only : c_int
use, intrinsic :: iso_fortran_env, only : int32
implicit none (type, external)

interface
integer(c_int) function cpu_count_c() bind(c)
import c_int
end function cpu_count_c
end interface

contains

module procedure get_cpu_count

get_cpu_count = int(cpu_count_c(), int32)

end procedure get_cpu_count


end submodule sysdet

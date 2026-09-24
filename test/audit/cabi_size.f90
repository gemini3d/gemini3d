function audit_params_size() bind(C) result(n)
use, intrinsic :: iso_c_binding, only: c_size_t,c_sizeof
use gemini3d, only: c_params
implicit none
type(c_params) :: p
integer(c_size_t) :: n
n=c_sizeof(p)
end function

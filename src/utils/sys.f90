module gemini3d_sysinfo
!! procedure to get information from the system

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, compiler_version

use filesystem, only : expanduser

implicit none (type, external)

private
public :: get_compiler_vendor, expand_envvar

contains


character(5) function get_compiler_vendor() result(vendor)
character(80) :: cvers
integer :: i, j
character(*), parameter :: vendors(2) = [character(5) :: "Intel", "GCC"]

cvers = compiler_version()

do j = 1,size(vendors)
  vendor = vendors(j)
  i = index(cvers, vendor)
  if (i > 0) exit
end do

if(vendor=="GCC") then
  vendor = "GNU"
elseif(i == 0) then
  vendor = ""
  write(stderr,'(A,/,A)') "could not determine compiler vendor from",cvers
end if
end function get_compiler_vendor


function expand_envvar(path) result(expanded)
!! replace @...@ string like metabuild system e.g. CMake, based on environment variable.
!!
!! NOTE: only expands the first @envvar@ substring. Nest calls if mutliple @envvar@ substrings

character(*), intent(in) :: path

character(:), allocatable :: expanded, envvar

integer :: i0, i1
integer :: L, istat
character(1000) :: buf

expanded = expanduser(path)

i0 = index(path, "@")
if (i0 < 1) return
i0 = i0

i1 = index(path(i0+1:), "@")
if (i1 < 1) return  !< a single @ without a matching @
i1 = i0 + i1

envvar = path(i0+1:i1-1)
if(len_trim(envvar) == 0) return  !< only blanks in envvar

call get_environment_variable(envvar, buf, length=L, status=istat)
if(istat /= 0) error stop "config:expand_envvar: environment variable not defined: " // envvar
if(L < 1) error stop "config:expand_envvar: environment variable empty: " // envvar

expanded = path(:i0-1) // trim(adjustl(buf)) // path(i1+1:)
end function expand_envvar

end module gemini3d_sysinfo

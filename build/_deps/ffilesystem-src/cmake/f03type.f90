module dummy
type :: path_t
private
character(:), allocatable :: path_str
contains
final :: destructor
end type

interface path_t
module procedure set_path
end interface

contains

subroutine destructor(self)
type(path_t), intent(inout) :: self
if(allocated(self%path_str)) deallocate(self%path_str)
end subroutine destructor

type(path_t) function set_path(path)
character(*), intent(in) :: path
allocate(character(100) :: set_path%path_str)
set_path%path_str = path
end function set_path

end module dummy


program main
use dummy
implicit none
type(path_t) :: p
p = path_t('.')
end program

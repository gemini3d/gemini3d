module nc4fortran
!! this is a dummy interface that errors intentionally

implicit none (type, external)

type netcdf_file
character(:), allocatable  :: filename
contains
procedure, public :: open, read, write, exist, close
end type netcdf_file

contains

subroutine open(self,filename, action,comp_lvl, verbose, debug)
class(netcdf_file), intent(inout)  :: self
character(*), intent(in)           :: filename
character(*), intent(in), optional :: action
integer, intent(in), optional      :: comp_lvl
logical, intent(in), optional :: verbose, debug
error stop 'NetCDF4 / nc4fortran not available'
end subroutine open


subroutine read(self, dname, value)
class(netcdf_file), intent(in)     :: self
character(*), intent(in)         :: dname
class(*), intent(inout)      :: value(..)
error stop 'NetCDF4 / nc4fortran not available'
end subroutine read

subroutine write(self, dname, value, dims)
class(netcdf_file), intent(in)     :: self
character(*), intent(in)         :: dname
class(*), intent(in)      :: value(..)
character(*), intent(in), optional :: dims(:)
error stop 'NetCDF4 / nc4fortran not available'
end subroutine

logical function exist(self, dname)
class(netcdf_file), intent(in)     :: self
character(*), intent(in)         :: dname
exist = .false.
error stop 'nc4fortran not available'
end function exist

subroutine close(self)
class(netcdf_file), intent(in) :: self
error stop 'NetCDF4 / nc4fortran not available'
end subroutine close

end module nc4fortran

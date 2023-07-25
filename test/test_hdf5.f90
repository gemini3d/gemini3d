program HDF5_standalone

use hdf5, only : HID_T, HSIZE_T, H5_INTEGER_KIND, h5kind_to_type, h5open_f, h5close_f, h5fclose_f, &
    h5fcreate_f, H5F_ACC_TRUNC_F, H5get_libversion_f
use h5lt, only : h5ltmake_dataset_f

implicit none

integer :: i, p
integer(HID_T) :: lid
character(*), parameter :: filename='test_minimal.h5'
integer :: major, minor, release

p = 42


!! check that repeated calls to h5open_f() do not cause problems as per docs
!! not calling h5open_f() at all makes failures as library isn't initialized
!! unlike C HDF5, Fortran HDF5 does not auto-initialize.
call h5open_f(i)
if(i /= 0) error stop "ERROR:hdf5_standalone_fortran: h5open_f failed, call #1"

call h5open_f(i)
if(i /= 0) error stop "ERROR:hdf5_standalone_fortran: h5open_f failed, call #2"

call h5open_f(i)
if(i /= 0) error stop "ERROR:hdf5_standalone_fortran: h5open_f failed, call #3"

call H5get_libversion_f(major, minor, release, i)
if (i /= 0) error stop "ERROR:hdf5_standalone_fortran:H5get_libversion: could not get HDF5 library version"

print '(a,i0,a1,i0,a1,i0)', "hdf5_standalone_fortran: HDF5 library version ", major, ".", minor, ".", release


call h5fcreate_f(filename, H5F_ACC_TRUNC_F, lid, i)
if (i/=0) error stop 'ERROR:hdf5_standalone_fortran: could not create file'
print *, 'hdf5_standalone_fortran: created '//filename

call h5ltmake_dataset_f(lid, "A", rank(p), shape(p, kind=HSIZE_T), h5kind_to_type(kind(p),H5_INTEGER_KIND), p, i)
if (i/=0) error stop 'ERROR:hdf5_standalone_fortran: could not create dataset A'
print *, 'hdf5_standalone_fortran: created variable'

call h5fclose_f(lid, i)
if (i/=0) error stop 'ERROR:hdf5_standalone_fortran: could not close file'
print *, 'hdf5_standalone_fortran: closed '//filename

call H5close_f(i)
if (i /= 0) error stop 'ERROR:hdf5_standalone_fortran: could not close hdf5 library'

! this is a Fortran-standard way to delete files
open(newunit=i, file=filename)
close(i, status='delete')

print *, 'OK: hdf5_standalone_fortran'

end program

!! This tests various text output formats.
!! It necessarily requires visual inspection, or run from an enclosing program with Regex,
!! such as a Python program configured to test as desired.
use, intrinsic:: iso_fortran_env, only: sp=>real32, dp=>real64
use formats, only: utsec2filename

implicit none

print *, 'Single precision lacks adequate precision for dates beyond millisecond'
print *, utsec2filename([2015,4,13], 12345.000003_sp)

if (utsec2filename([2015,4,13], 12345.678910_dp) /= '20150413_12345.678910.dat') error stop 'mismatch real64'

if (utsec2filename([2015,4,13],  345.678910_dp) /= '20150413_00345.678910.dat') error stop 'mismatch leading zeros'

end program
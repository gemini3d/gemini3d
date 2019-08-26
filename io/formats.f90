module date_formats

use, intrinsic :: iso_fortran_env, only: real64, real32, int32, int64

implicit none

contains

pure function date_filename(outdir,ymd,UTsec)  result(filename)
!! GENERATE A FILENAME STRING OUT OF A GIVEN DATE/TIME

character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
character(:), allocatable :: filename


!> assemble
filename = outdir // '/' // utsec2filename(ymd, UTsec)

end function date_filename


pure character(25) function utsec2filename(ymd, UTsec) result(fn)
!! file name is exactly 25 characters long, per Matt Z's de facto spec.
!! FIXME: until we go to integer UTsec (microsec) we round to nearest microsecond
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
!! UTC second: real [0.0 .. 86400.0)
integer(int64) :: usec, seconds, frac
character(12) :: sec_str

select type(UTsec)
  type is (real(real32))
    usec = nint(UTsec * 1000000_int64, int64)
  type is (real(real64))
    usec = nint(UTsec * 1000000_int64, int64)
  type is (integer(int32))
    usec = int(UTsec, int64) * 1000000_int64
  type is (integer(int64))
    usec = UTsec * 1000000_int64
  class default
    error stop "io/formats.f90:utsec2filename unknown UTsec type"
end select

seconds = usec / 1000000_int64
if (seconds < 0 .or. seconds >= 86400) error stop 'io/formats.f90:utsec2filename did NOT satisfy 0 <= seconds < 86400'
frac = usec-seconds*1000000

write(sec_str, '(I5.5, A1, I6.6)') seconds, '.', frac

write(fn,'(i4,2I2.2,a17)') ymd(1), ymd(2:), '_' // sec_str // '.dat'

end function utsec2filename

end module date_formats

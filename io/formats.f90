module formats

use, intrinsic :: iso_fortran_env, only: real64, real32, int64

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


pure character(25) function utsec2filename(ymd, utsec) result(fn)
!! file name is exactly 25 characters long, per Matt Z's de facto spec.
!! FIXME: until we go to integer UTsec (microsec) we round to nearest microsecond
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
!! UTC second: real [0.0 .. 86400.0)
integer(int64) :: usec, seconds, frac
character(12) :: sec_str

select type(UTsec)
  type is (real(real32))
    usec = nint(UTsec * 1000000, int64)
  type is (real(real64))
    usec = nint(UTsec * 1000000, int64)
  class default
    error stop "unknown type into utsec2filename()"
end select

seconds = usec / 1000000
frac = usec-seconds*1000000

if (seconds < 0 .or. seconds >= 86400) error stop '0 <= seconds < 86400'
write(sec_str, '(I5.5, A1, I6.6)') seconds, '.', frac

write(fn,'(i4,2I2.2,a17)') ymd(1), ymd(2:), '_' // sec_str // '.dat'

end function utsec2filename

end module formats

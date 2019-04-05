module formats

use, intrinsic :: iso_fortran_env, only: real64, real32

implicit none

contains

pure function date_filename(outdir,ymd,UTsec)
!! GENERATE A FILENAME STRING OUT OF A GIVEN DATE/TIME

character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
character(:), allocatable :: date_filename


!> assemble
date_filename = outdir // '/' // utsec2filename(ymd, UTsec)

end function date_filename


pure character(25) function utsec2filename(ymd, utsec) result(fn)
!! file name is exactly 25 characters long, per Matt Z's de facto spec.
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
!! UTC second: real [0.0 .. 86400.0)

character(12) :: sec_str

select type(UTsec)
  type is (real(real32))
    write(sec_str, '(I5.5, F7.6)') int(UTsec), UTsec-int(UTsec)
  type is (real(real64))
    write(sec_str, '(I5.5, F7.6)') int(UTsec), UTsec-int(UTsec)
end select

write(fn,'(i4,2I2.2,a17)') ymd, '_' // sec_str // '.dat'

end function utsec2filename

end module formats
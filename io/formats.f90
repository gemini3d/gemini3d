module formats

use, intrinsic :: iso_fortran_env, only: real64, real32

implicit none

interface date_filename
  module procedure old_date_filename
end interface date_filename

contains

pure function new_date_filename(outdir,ymd,UTsec)
!! GENERATE A FILENAME STRING OUT OF A GIVEN DATE/TIME

character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
character(:), allocatable :: new_date_filename


!> assemble
new_date_filename = outdir // '/' // utsec2filename(ymd, UTsec)

end function new_date_filename


pure character(25) function utsec2filename(ymd, utsec) result(fn)
!! file name is exactly 25 characters long, per Matt Z's de facto spec.
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
!! UTC second: real [0.0 .. 86400.0)

!! the rather verbose handling below is to ensure proper formatting across all compilers and OS.

character(12) :: sec_str

select type(UTsec)
  type is (real(real32))
    if (UTsec < 0 .or. UTsec >= 86400) error stop '0 <= utsec < 86400'
    write(sec_str, '(I5.5, A1, I6.6)') int(UTsec), '.', int((UTsec-int(UTsec))*1000000)
  type is (real(real64))
    if (UTsec < 0 .or. UTsec >= 86400) error stop '0 <= utsec < 86400'
    write(sec_str, '(I5.5, A1, I6.6)') int(UTsec), '.', int((UTsec-int(UTsec))*1000000)
end select

write(fn,'(i4,2I2.2,a17)') ymd(1), ymd(2:), '_' // sec_str // '.dat'

end function utsec2filename


!OLD LEGACY CODE FOR MATT Z.
  function old_date_filename(outdir,ymd,UTsec)

    !------------------------------------------------------------
    !-------GENERATE A FILENAME STRING OUT OF A GIVEN DATE/TIME
    !------------------------------------------------------------

    character(*), intent(in) :: outdir
    integer, dimension(3), intent(in) :: ymd
    real(8), intent(in) :: UTsec
    character(:), allocatable :: old_date_filename

    integer :: ldigits,idigits
    character(256) :: filename,tmpchar,tmpchar2
    character(512) :: tmpfilename,filenamefull


    !FORM OUTPUT FILENAME BASED ON DATE AND TIME (this was unbelievably squirrely to work out)
    write(filename,'(f12.6,a4)') UTsec,'.dat'    !file name that has 6 decimal points on time stamp
    filename=adjustl(filename)                   !slam the chars. to the left and remove trailing blanks
    ldigits=5                                    !pad the filename with the appropriate number zero characters
    if (UTsec<1d0) then
      idigits=1
    else
      idigits=floor(log10(UTsec))+1
    end if
    tmpchar=filename
    do while(idigits<ldigits)
      write(tmpchar2,*) '0',trim(tmpchar)
      tmpchar=adjustl(tmpchar2)      
      idigits=idigits+1
    end do
    filename=tmpchar

    !day
    write(tmpchar,*) ymd(3)
    tmpchar=adjustl(tmpchar)
    if (ymd(3)<10) then
      write(tmpchar2,*) '0',trim(tmpchar)
      tmpchar=adjustl(tmpchar2)
    end if
    !write is dumb and doesn't recognize previous trims...  I hate string manipulation...
    write(tmpfilename,*) trim(tmpchar),'_',trim(filename)
    tmpfilename=adjustl(tmpfilename)
    filename=tmpfilename(1:256)

    !month
    write(tmpchar,*) ymd(2)
    tmpchar=adjustl(tmpchar)
    if (ymd(2)<10) then
      write(tmpchar2,*) '0',trim(tmpchar)
      tmpchar=adjustl(tmpchar2)
    end if
    write(tmpfilename,*) trim(tmpchar),trim(filename)
    tmpfilename=adjustl(tmpfilename)
    filename=tmpfilename(1:256)

    !year
    write(tmpchar,*) ymd(1)
    tmpchar=adjustl(tmpchar)
    write(tmpfilename,*) trim(tmpchar),trim(filename)
    tmpfilename=adjustl(tmpfilename)
    filename=tmpfilename(1:256)
    write(filenamefull,*) outdir,'/',trim(filename)
    old_date_filename=trim(adjustl(filenamefull))

  end function old_date_filename


end module formats

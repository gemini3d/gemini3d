program msis_driver
!! will write to stdout if "-" specified, so we avoid printing to console unless file output is used
use, intrinsic:: iso_fortran_env, only: sp=>real32, stderr=>error_unit, stdout=>output_unit, stdin=>input_unit

implicit none (type, external)

integer, parameter :: mass=48
integer :: iyd,sec,lz, i
real(sp) :: f107a,f107,ap(7),stl,apday,ap3
real(sp) :: d(9),t(2)
real(sp), allocatable :: glat(:),glon(:),alt(:)
integer :: u
character(256) :: buf
character(:), allocatable :: infile,outfile

external :: meters, gtd7

!> read in msis inputs
if (command_argument_count() < 2) error stop 'msis_setup: must specify input and output filenames'

call get_command_argument(1,buf)
infile = trim(buf)

if (infile == "-" .or. infile(len(infile)-3:) == ".txt") then
  call get_text_input(infile, iyd,sec,f107a,f107,apday,ap3,lz, glat, glon, alt)
else
  call get_binary_input(infile,iyd,sec,f107a,f107,apday,ap3,lz, glat, glon, alt)
endif

!> Run MSIS
ap(1:7)=apday
ap(2)=ap3

!> switch to mksa units
call meters(.true.)

!> output file
call get_command_argument(2, buf)
outfile = trim(buf)

if (outfile == '-') then
  u = stdout
elseif(outfile(len(outfile)-3:) == ".txt") then
  open(newunit=u, file=outfile, status='replace',action='write')
else
  open(newunit=u,file=outfile,status='replace',form='unformatted',access='stream', action='write')
endif

!> call to msis routine
do i=1,lz
  stl = sec/3600. + glon(i)/15.
  call gtd7(iyd, real(sec, sp),alt(i),glat(i),glon(i),stl,f107a,f107,ap,mass,d,t)

  if (outfile == '-') then
    write(u,'(F9.2, 9ES15.6, F9.2)') alt(i),d(1:9),t(2)
  elseif(outfile(len(outfile)-3:) == ".txt") then
    write(u,'(F9.2, 9ES15.6, F9.2)') alt(i),d(1:9),t(2)
  else
    write(u) alt(i),d(1:9),t(2)
  endif
end do

if (outfile == '-') stop

close(u)

inquire(file=outfile, size=i)
print *,'msis_setup: wrote ',i,' bytes to ',outfile
if (i==0) error stop 'msis_setup failed to write file'

contains

subroutine get_text_input(filename, iyd,sec,f107a,f107,apday,ap3,lz, glat, glon, alt)

character(*), intent(in) :: filename
integer, intent(out) :: iyd,sec,lz
real(sp), intent(out) :: f107a,f107,apday,ap3
real(sp), intent(out), allocatable :: glat(:),glon(:),alt(:)
integer :: u, i

if(filename == "-") then
  u = stdin
else
  open(newunit=u, file=filename, status="old", action="read")
endif

read(u, *, iostat=i) iyd
if (i/=0) error stop "iyd: integer yyddd"
read(u, *, iostat=i) sec
if (i/=0) error stop "sec: integer seconds since midniight"
read(u, *, iostat=i) f107a, f107, apday, ap3
if (i/=0) error stop "expecting: f107a, f107, apday, ap3"
read(u, *, iostat=i) lz
if (i/=0) error stop "lz: expecting integer number of altitudes"

call check_lz(lz)
allocate(glat(lz),glon(lz),alt(lz))
read(u,*, iostat=i) glat
if (i/=0)then
  write(stderr,*) 'read iostat', i
  error stop "glat: ran out of input elements unexpectedly"
endif

read(u,*, iostat=i) glon
if (i/=0) then
  write(stderr,*) 'read iostat', i
  error stop "glon: ran out of input elements unexpectedly"
endif
read(u,*, iostat=i) alt
if (i/=0) error stop "alt: ran out of input elements unexpectedly"

if (filename /= "-") close(u)

end subroutine get_text_input


subroutine get_binary_input(filename,iyd,sec,f107a,f107,apday,ap3,lz, glat, glon, alt)

character(*), intent(in) :: filename
integer, intent(out) :: iyd,sec,lz
real(sp), intent(out) :: f107a,f107,apday,ap3
real(sp), intent(out), allocatable :: glat(:),glon(:),alt(:)

integer :: u
open(newunit=u,file=infile, status='old',form='unformatted',access='stream', action='read')
!! use binary to reduce file size and read times

read(u) iyd
read(u) sec
read(u) f107a
read(u) f107
read(u) apday
read(u) ap3
read(u) lz

call check_lz(lz)
allocate(glat(lz),glon(lz),alt(lz))
read(u) glat,glon,alt

close(u)

end subroutine get_binary_input


subroutine check_lz(lz)

integer, intent(in) :: lz
character(256) :: buf
integer :: i

if (lz<1) error stop 'lz must be positive'

call get_command_argument(3, buf, status=i)
if (i==0) then
  read(buf,*) i
  if (i /= lz) then
    write(stderr,*) 'expected ',i,' grid points but read ',lz
    error stop
  endif
endif

end subroutine check_lz


end program

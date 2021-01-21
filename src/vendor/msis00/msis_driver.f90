program msis_driver
!! will write to stdout if "-" specified, so we avoid printing to console unless file output is used

use phys_consts, only : comp_lvl
use msis_interface, only : msis_gtd7, msis_gtd8, msisinit
use h5fortran, only : hdf5_file, hsize_t
use, intrinsic:: iso_fortran_env, only : real32, stderr=>error_unit, stdout=>output_unit, stdin=>input_unit

implicit none (type, external)

integer :: i,j,k
real(real32) :: doy,sec,f107a,f107, Ap(7)
real(real32), allocatable, dimension(:,:,:,:) :: Dn, Tn
real(real32), allocatable, dimension(:,:,:) :: glat, glon, alt
integer :: u, msis_version, lx1, lx2, lx3
character(256) :: buf
character(:), allocatable :: infile,outfile

!> user options
if (command_argument_count() < 2) error stop 'msis_setup: must specify input and output filenames'

call get_command_argument(1,buf)
infile = trim(buf)
call get_command_argument(2, buf)
outfile = trim(buf)

msis_version = 0
if(command_argument_count() > 2) then
  call get_command_argument(3, buf)
  read(buf, '(I3)') msis_version
endif

!> select input format
call input_hdf5(infile,doy,sec,f107a,f107,Ap, glat, glon, alt)

!> Run MSIS
lx1 = size(alt,1)
lx2 = size(alt,2)
lx3 = size(alt,3)
allocate(Dn(lx1,lx2,lx3, 9), Tn(lx1,lx2,lx3, 2))

if(msis_version == 20) call msisinit(parmfile='msis20.parm')

do i=1,lx1
  do j=1,lx2
    do k=1,lx3

      if(msis_version == 0) then
        call msis_gtd7(doy, sec, alt(i,j,k), glat(i,j,k), glon(i,j,k), f107a, f107, Ap, Dn(i,j,k,:), Tn(i,j,k,:), use_meters=.true.)
      elseif(msis_version == 20) then
        call msis_gtd8(doy, sec, alt(i,j,k), glat(i,j,k), glon(i,j,k), f107a, f107, Ap, Dn(i,j,k,:), Tn(i,j,k,:))
      else
        error stop 'expected msis_version = {0,20}'
      end if

    end do
  end do
end do

call output_hdf5(outfile, alt, Dn, Tn)

contains


subroutine input_hdf5(filename,doy,sec,f107a,f107,Ap7, glat, glon, alt)
!! use binary to reduce file size and read times
character(*), intent(in) :: filename
real(real32), intent(out) :: doy,sec,f107a,f107,ap7(7)
real(real32), intent(out), allocatable :: glat(:,:,:), glon(:,:,:), alt(:,:,:)

type(hdf5_file) :: hf
integer(hsize_t), allocatable :: dims(:)
integer:: lx1,lx2,lx3

call hf%initialize(filename, status='old', action='read')

call hf%read("/doy", doy)
call hf%read("/UTsec", sec)
call hf%read("/f107a", f107a)
call hf%read("/f107", f107)
call hf%read("/Ap", Ap7)

call hf%shape("/glat", dims)
lx1 = int(dims(1))
lx2 = int(dims(2))
lx3 = int(dims(3))

allocate(glat(lx1,lx2,lx3), glon(lx1,lx2,lx3), alt(lx1,lx2,lx3))

call hf%read("/glat", glat)
call hf%read("/glon", glon)
call hf%read("/alt", alt)

call hf%finalize()

end subroutine input_hdf5


subroutine output_hdf5(filename, alt, Dn, Tn)

character(*), intent(in) :: filename
real(real32), intent(in) :: alt(:,:,:), Dn(:,:,:,:), Tn(:,:,:,:)
type(hdf5_file) :: hf

call hf%initialize(filename, status="replace", action="write",comp_lvl=comp_lvl)

call hf%write("/alt", alt)

call hf%write("/nHe", Dn(:,:,:,1))
call hf%write("/nO", Dn(:,:,:,2))
call hf%write("/nN2", Dn(:,:,:,3))
call hf%write("/nO2", Dn(:,:,:,4))
call hf%write("/nAr", Dn(:,:,:,5))
call hf%write("/TotalMassDensity", Dn(:,:,:,6))
call hf%write("/nH", Dn(:,:,:,7))
call hf%write("/nN", Dn(:,:,:,8))
call hf%write("/nOana", Dn(:,:,:,9))

call hf%write("/Tn", Tn(:,:,:,2))
call hf%write("/Texo", Tn(:,:,:,1))


call hf%finalize()

end subroutine output_hdf5


! ------- below not used anymore

subroutine output_text(filename, alt, Dn, Tn)

character(*), intent(in) :: filename
real(real32), intent(in) :: alt(:), Dn(:,:), Tn(:,:)
logical :: pipe

pipe = filename == '-'

if (pipe) then
  u = stdout
else
  open(newunit=u, file=outfile, status='replace',action='write')
endif

do i = 1,size(alt)
  write(u,'(F9.2, 9ES15.6, F9.2)') alt(i), Dn(i,:), Tn(i,2)
end do

if (outfile /= '-') close(u)

end subroutine output_text


subroutine input_text(filename, doy,sec,f107a,f107,apday,ap3, glat, glon, alt)

character(*), intent(in) :: filename
real(real32), intent(out) :: doy,sec,f107a,f107,apday,ap3
real(real32), intent(out), allocatable :: glat(:),glon(:),alt(:)

integer :: u, i, lz

if(filename == "-") then
  u = stdin
else
  open(newunit=u, file=filename, status="old", action="read")
endif

read(u, *, iostat=i) doy
if (i/=0) error stop "doy: day of year (1..366)"
read(u, *, iostat=i) sec
if (i/=0) error stop "sec: seconds since UTC midnight"
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

end subroutine input_text


subroutine check_lz(lz)

integer, intent(in) :: lz
character(256) :: buf
integer :: i

if (lz < 1) error stop 'lz must be positive'

call get_command_argument(4, buf, status=i)
if (i==0) then
  read(buf, *) i
  if (i /= lz) then
    write(stderr,*) 'expected ',i,' grid points but read ',lz
    error stop
  endif
endif

end subroutine check_lz


end program

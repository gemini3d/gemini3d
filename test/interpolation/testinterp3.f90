program testinterp3

use phys_consts, only: wp,pi
use interpolation, only : interp3
use h5fortran, only: hdf5_file

implicit none (type, external)

character(1024) :: argv
type(hdf5_file) :: hout
integer, parameter :: lx1=80, lx2=90, lx3=100
integer, parameter :: lx1i=16, lx2i=96, lx3i=32
real, parameter :: stride=0.5

real(wp), allocatable, dimension(:) :: x1, x2, x3, x1i, x2i, x3i
real(wp), allocatable, dimension(:,:,:) :: f, fi

real(wp), allocatable, dimension(:) :: x1ilist,x2ilist,x3ilist,filist

integer :: ix1,ix2,ix3,ik, ierr


!> grid for original data
x1=[ ((real(ix1,wp)-1)*stride, ix1=1,lx1) ]
x2=[ ((real(ix2,wp)-1)*stride, ix2=1,lx2) ]
x3=[ ((real(ix3,wp)-1)*stride, ix3=1,lx3) ]

!> center grid points at zero
x1=x1-sum(x1)/size(x1,1)
x2=x2-sum(x2)/size(x2,1)
x3=x3-sum(x3)/size(x3,1)

allocate(f(lx1,lx2,lx3))
!> test function
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      f(ix1,ix2,ix3)=sin(2*pi/5._wp*x1(ix1))*cos(2*pi/20._wp*x2(ix2))*sin(2*pi/15._wp*x3(ix3))
    end do
  end do
end do

!> grid for interpolated data
x1i = [ ((real(ix1)-1) * stride/(lx1i/real(lx1)), ix1=1,lx1i) ]
x2i = [ ((real(ix2)-1) * stride/(lx2i/real(lx2)), ix2=1,lx2i) ]
x3i = [ ((real(ix3)-1) * stride/(lx3i/real(lx3)), ix3=1,lx3i) ]

!> center grid points at zero
x1i=x1i-sum(x1i)/size(x1i,1)
x2i=x2i-sum(x2i)/size(x2i,1)
x3i=x3i-sum(x3i)/size(x3i,1)

allocate(x1ilist(lx1i*lx2i*lx3i))
allocate(x2ilist, x3ilist, filist, mold=x1ilist)

!> try a 3d interpolation
do ix3=1,lx3i
  do ix2=1,lx2i
    do ix1=1,lx1i
      ik=(ix3-1)*lx1i*lx2i+(ix2-1)*lx1i+ix1
      x1ilist(ik)=x1i(ix1)
      x2ilist(ik)=x2i(ix2)
      x3ilist(ik)=x3i(ix3)
    end do
  end do
end do

print *, "interp3d: Starting test interpolation"
print *, shape(x1), shape(x2), shape(x3), shape(f)
print *, shape(x1ilist), shape(x2ilist), shape(x3ilist), shape(filist)
filist = interp3(x1, x2, x3, f, x1ilist, x2ilist, x3ilist)

fi = reshape(filist, [lx1i,lx2i,lx3i])

!> sanity check
if (any(shape(fi) /= [lx1i, lx2i, lx3i])) error stop "test_interp3d: not expected shape"


call get_command_argument(1, argv, status=ierr)
if(ierr /= 0) error stop 'please specify input filename'


print "(A,/,A,/,A)", "interp3d: Finished test interpolation"
!> dump results to a file so we can check things
call hout%open(trim(argv), action="w")

call hout%write("/lx1", lx1)
call hout%write("/lx2", lx2)
call hout%write("/lx3", lx3)
call hout%write("/x1", x1)
call hout%write("/x2", x2)
call hout%write("/x3", x3)
call hout%write("/f", f)

call hout%close()


call get_command_argument(2, argv)
if(argv=="") error stop 'please specify output filename'

call hout%open(trim(argv), action="w")

call hout%write("/lx1", lx1i)
call hout%write("/lx2", lx2i)
call hout%write("/lx3", lx3i)
call hout%write("/x1", x1i)
call hout%write("/x2", x2i)
call hout%write("/x3", x3i)
call hout%write("/f", fi)

call hout%close()

end program


!> has no problem with > 2GB output files
! block
!   integer :: u
!   open(newunit=u,file='input3D.dat',status='replace',form='unformatted',access='stream', action='write')
!   write(u) lx1,lx2,lx3
!   write(u) x1,x2,x3,f
!   close(u)
! end block

! block
!   integer :: u
!   open(newunit=u,file='output3D.dat',status='replace',form='unformatted',access='stream', action='write')
!   write(u) lx1i,lx2i,lx3i
!   write(u) x1i,x2i,x3i,fi   !< since only interpolating in x1
!   close(u)
! end block

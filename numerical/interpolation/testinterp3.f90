program testinterp3

use phys_consts, only: wp,pi
use interpolation
use h5fortran, only: hdf5_file

implicit none

type(hdf5_file) :: hout
integer, parameter :: lx1=80, lx2=90, lx3=100
integer, parameter :: lx1i=256, lx2i=256, lx3i=256
real(wp), parameter :: stride=0.5_wp

real(wp) :: x1(lx1),x2(lx2),x3(lx3),f(lx1,lx2,lx3)
real(wp) :: x1i(lx1i),x2i(lx2i),x3i(lx3i),fi(lx1i,lx2i,lx3i)

real(wp), dimension(1:lx1i*lx2i*lx3i) :: x1ilist,x2ilist,x3ilist,filist

integer :: ix1,ix2,ix3,ik, ierr


!grid for original data
x1=[ ((real(ix1,wp)-1._wp)*stride, ix1=1,lx1) ]
x2=[ ((real(ix2,wp)-1._wp)*stride, ix2=1,lx2) ]
x3=[ ((real(ix3,wp)-1._wp)*stride, ix3=1,lx3) ]

!center grid points at zero
x1=x1-sum(x1)/size(x1,1)
x2=x2-sum(x2)/size(x2,1)
x3=x3-sum(x3)/size(x3,1)

!test function
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      f(ix1,ix2,ix3)=sin(2._wp*pi/5._wp*x1(ix1))*cos(2._wp*pi/20._wp*x2(ix2))*sin(2._wp*pi/15._wp*x3(ix3))
    end do
  end do
end do

!> grid for interpolated data
x1i=[ ((real(ix1,wp)-1)*stride/(lx1i/lx1), ix1=1,lx1i) ]
x2i=[ ((real(ix2,wp)-1)*stride/(lx2i/lx2), ix2=1,lx2i) ]
x3i=[ ((real(ix3,wp)-1)*stride/(lx3i/lx3), ix3=1,lx3i) ]

!> center grid points at zero
x1i=x1i-sum(x1i)/size(x1i,1)
x2i=x2i-sum(x2i)/size(x2i,1)
x3i=x3i-sum(x3i)/size(x3i,1)


!> try a 333d interpolation
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
filist=interp3(x1,x2,x3,f,x1ilist,x2ilist,x3ilist)
fi=reshape(filist,[lx1i,lx2i,lx3i])


!> dump results to a file so we can check things
call hout%initialize("input3d.h5", ierr, status="replace", action="write")
if(ierr/=0) error stop 'interp3: could not open input file'

call hout%write("/lx1", lx1, ierr)
call hout%write("/lx2", lx2, ierr)
call hout%write("/lx3", lx3, ierr)
call hout%write("/x1", x1, ierr)
call hout%write("/x2", x2, ierr)
call hout%write("/x3", x3, ierr)
call hout%write("/f", f, ierr)

call hout%finalize(ierr)
if(ierr/=0) error stop 'interp3: could not close input file'


call hout%initialize("output3d.h5", ierr, status="replace", action="write")
if(ierr/=0) error stop 'interp3: could not open output file'

call hout%write("/lx1", lx1i, ierr)
call hout%write("/lx2", lx2i, ierr)
call hout%write("/lx3", lx3i, ierr)
call hout%write("/x1", x1i, ierr)
call hout%write("/x2", x2i, ierr)
call hout%write("/x3", x3i, ierr)
call hout%write("/f", fi, ierr)

call hout%finalize(ierr)
if(ierr/=0) error stop 'interp3: could not close output file'

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
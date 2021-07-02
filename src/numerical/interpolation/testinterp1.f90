program test_interp1
!! Need program statement for FORD
use phys_consts, only: wp,pi
use interpolation

implicit none (type, external)

integer, parameter :: lx1=25, lx2=50
integer, parameter :: lx1i=200, lx2i=400
real(wp), parameter :: stride=0.5_wp

integer :: ix1,ix2
integer :: u

real(wp), allocatable, dimension(:) :: x1, x2, x1i, x2i
real(wp), allocatable, dimension(:,:) :: f, fi  !< since only doing x1 interpolation

allocate(x1(lx1), x2(lx2), f(lx1,lx2), x1i(lx1i), x2i(lx2i), fi(lx1i,lx2i))

!grid for original data
x1=[ ((real(ix1,wp)-1._wp)*stride, ix1=1,lx1) ]
x2=[ ((real(ix2,wp)-1._wp)*stride, ix2=1,lx2) ]


!test function
do ix2=1,lx2
  f(:,ix2)=sin(2._wp*pi/5._wp*x1)*cos(2._wp*pi/5._wp*x2(ix2))
end do


!grid for interpolated data
x1i=[ ((real(ix1,wp)-1._wp)*stride/(lx1i/lx1), ix1=1,lx1i) ]
x2i=[ ((real(ix2,wp)-1._wp)*stride/(lx2i/lx2), ix2=1,lx2i) ]


!try a 1d interpolation for each x2 value used in simulation
do ix2=1,lx2
  fi(:,ix2)=interp1(x1,f(:,ix2),x1i)
end do

!dump results to a file so we can check things
open(newunit=u,file='input1D.dat',status='replace',form='unformatted',access='stream', action='write')
write(u) lx1,lx2
write(u) x1,x2,f
close(u)

open(newunit=u,file='output1D.dat',status='replace',form='unformatted',access='stream', action='write')
write(u) lx1i,lx2
write(u) x1i,x2,fi   !since only interpolating in x1
close(u)

end program

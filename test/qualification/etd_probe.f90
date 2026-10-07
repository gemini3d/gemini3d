! SPDX-License-Identifier: Apache-2.0
program etd_probe
use phys_consts, only: wp
use calculus, only: ETD_uncoupled
implicit none(type,external)
real(wp) :: f(1,1,1),p(1,1,1),l(1,1,1),y(1,1,1),dt
integer :: ios
do
  read(*,*,iostat=ios) f,p,l,dt
  if(ios<0) exit
  if(ios/=0) error stop 'Expected f P L dt'
  y=ETD_uncoupled(f,p,l,dt)
  write(*,'(ES26.17E3)') y
enddo
end program

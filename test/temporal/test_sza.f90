program test_sza
!! test day rollover
use phys_consts, only: wp
use timeutils, only: sza

implicit none (type, external)

real(wp), parameter :: pi = 4._wp*atan(1._wp)

print *, sza(2015,1,1, 0._wp, 0._wp, 0._wp)*180/pi
print *, sza(2015,6,21, 0._wp, 0._wp, 0._wp)*180/pi
print *, sza(2015,6,1, 43200._wp, 0._wp, 0._wp)*180/pi

print *, sza(2015,1,1, 0._wp, 45._wp, 45._wp)*180/pi
print *, sza(2015,6,21, 0._wp, 45._wp, 45._wp)*180/pi
print *, sza(2015,6,21, 43200._wp, 45._wp, 45._wp)*180/pi

end program

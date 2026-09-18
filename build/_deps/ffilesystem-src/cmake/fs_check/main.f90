program main

use filesystem_test

implicit none

integer :: i

i = merge(1, 0, has_filename("./"))

stop i

end program

program test_is_abs

use filesystem

logical :: a

if (is_absolute("")) error stop "blank is not absolute"

a = is_absolute("/")

if (is_windows()) then
  if (a) error stop "is_absolute(/) on Windows should be false"
else
  if (.not. a) error stop "is_absolute(/) on Unix should be true"
end if

end program

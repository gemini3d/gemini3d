submodule (pathlib) pathlib_gcc

implicit none (type, external)

contains

module procedure directory_exists
!! For GCC Gfortran, similar for other compilers
integer :: i, statb(13)
character(:), allocatable :: wk

wk = expanduser(path)

!! must not have trailing slash on Windows
i = len_trim(wk)
if (wk(i:i) == char(92) .or. wk(i:i) == '/') wk = wk(1:i-1)


inquire(file=wk, exist=exists)
if(.not.exists) return

i = stat(wk, statb)
if(i /= 0) then
  exists = .false.
  return
endif

i = iand(statb(3), O'0040000')
exists = i == 16384

! print '(O8)', statb(3)

end procedure directory_exists

end submodule pathlib_gcc

submodule (pathlib) pathlib_gcc

implicit none (type, external)

contains

module procedure is_dir

integer :: i, statb(13)
character(:), allocatable :: wk

wk = expanduser(path)

!! must not have trailing slash on Windows
i = len_trim(wk)
if (wk(i:i) == '/') wk = wk(1:i-1)

inquire(file=wk, exist=is_dir)
if(.not.is_dir) return

i = stat(wk, statb)
if(i /= 0) then
  is_dir = .false.
  return
endif

i = iand(statb(3), O'0040000')
is_dir = i == 16384

! print '(O8)', statb(3)

end procedure is_dir

end submodule pathlib_gcc

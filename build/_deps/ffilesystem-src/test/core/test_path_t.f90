!> use a dummy module instead of
!> program-contains to workaround bug in NVHPC through at least 24.3

module dummy
use filesystem, only : path_t, devnull
implicit none

contains

subroutine test_setter_getter()

type(path_t) :: p1

p1 = path_t("a/b/c")

if (p1%path(2,3) /= "/b") error stop "getter start,end"
if (p1%path(3,3) /= "b") error stop "getter same"
if (p1%path(2) /= "/b/c") error stop "getter start only"

end subroutine

end module


program path_methods

use, intrinsic :: iso_fortran_env, only : stderr => ERROR_UNIT, int64
use filesystem
use dummy, only : test_setter_getter

implicit none

call test_setter_getter()
print '(a)', "OK: setter/getter"

valgrind : block

type(path_t) :: p1, p2
integer :: L1
integer(int64) :: i64
character(:), allocatable :: s1, s2
character(9) :: c9
logical :: b

!> as_posix
if(is_windows()) then
  s2 = "C:" // char(92) // "a" // char(92) // "b"
  p1 = path_t(s2)
  p1 = p1%as_posix()
  s1 = "C:/a/b"
else
  s2 = "/a/b"
  p1 = path_t(s2)
  p1 = p1%as_posix()
  s1 = s2
end if
if(p1%path() /= s1) then
  write(stderr, '(a,i0,a,i0)') "as_posix(" // s2 // "): expected " // s1 // " length: ", len(s1), &
    " but got " // p1%path() // " length:", p1%length()
  error stop
endif

print '(a)', "PASS: %as_posix(" // s2 // ") = " // p1%path()

!> absolute
p1 = path_t(s2)
p2 = path_t("a/b")
if(.not. p1%is_absolute()) error stop "is_absolute(abs)"
if(p2%is_absolute()) error stop "is_absolute(rel)"

!> expanduser
p1 = path_t("~")
p1 = p1%expanduser()
if (p1%path() /= get_homedir()) error stop "expanduser"

!> canonical
p1 = path_t(".")
p1 = p1%canonical()
L1 = p1%length()
if (L1 == 0) then
  write(stderr,*) "ERROR: canonical(.) " // p1%path()
  error stop
end if

if (.not. same_file(p1%path(), get_cwd())) then
  write(stderr,*) "ERROR: canonical(.) " // p1%path() // " /= get_cwd: " // get_cwd()
  error stop
end if

!> resolve
p1 = path_t(".")
p1 = p1%resolve()
L1 = p1%length()
if (L1 == 0) then
  write(stderr,*) "ERROR: resolve(.) " // p1%path()
  error stop
end if

if (.not. same_file(p1%path(), get_cwd())) then
  write(stderr,*) "ERROR: resolve(.) " // p1%path() // " /= get_cwd: " // get_cwd()
  error stop
end if

!> filename
p1 = path_t("a/b/c")
if (p1%file_name() /= "c") error stop "filename"

!> is_dir, is_file, exists
s2 = getarg(0)
if (len_trim(s2) == 0) error stop "could not get own program name"

p1 = path_t(s2)
if (.not. p1%exists()) error stop "ERROR:test: exists(file)"
if (p1%is_dir()) error stop "ERROR:test: is_dir(file)"
if (.not. p1%is_file()) error stop "ERROR:test: is_file(file)"
p1 = path_t(p1%parent())
if (.not. p1%is_dir()) error stop "ERROR:test: is_dir(dir)"

!> is_exe
p1 = path_t(s2)
if (.not. p1%is_exe()) error stop "ERROR:test: is_exe(exe)"

!> file_size
p1 = path_t(s2)
if (p1%file_size() == 0) error stop "ERROR:test: file_size"

!> modtime
i64 = p1%modtime()
if (i64 == 0) error stop "ERROR:test: modtime"

!> get_permissions
c9 = p1%get_permissions()
print '(a)', "get_permissions: " // c9
if(c9(1:3) /= "rwx") error stop "ERROR:test: get_permissions expected user rwx"

!> is_subdir
p1 = path_t("a/b")
if (.not.p1%is_subdir("a")) error stop "ERROR:test: is_subdir(a)"
if (p1%is_subdir("a/b/c")) error stop "ERROR:test: is_subdir(a/b/c)"

!> join
p1 = path_t("a")
p1 = p1%join("b")
if (p1%path() /= "a/b") error stop "ERROR:test: join"

!> normal
p1 = path_t("a/b/../c")
p1 = p1%normal()
if (p1%path() /= "a/c") error stop "ERROR:test: normal"

!> parent
p1 = path_t("a/b/c")
if (p1%parent() /= "a/b") error stop "ERROR:test: parent"

!> root
p1 = path_t("/c/b")
s1 = "/"
if (p1%root() /= s1) error stop "ERROR:test: root: " // p1%root() // " /= " // s1

!> stem
p1 = path_t("a/b.txt")
if (p1%stem() /= "b") error stop "ERROR:test: stem: " // p1%stem()

!> suffix
s1 = p1%suffix()
if (s1 /= ".txt") then
  write(stderr,'(a)') "ERROR:test: suffix: " // s1 // " /= .txt"
  error stop
end if
print '(a)', "PASS: %suffix: " // s1

!> touch
p1 = path_t("touch.txt")
call p1%touch()
if (.not. p1%is_file()) error stop "touch"

!> copy_file, with_suffix
p2 = p1%with_suffix(".dat")
if(p2%path() /= "touch.dat") error stop "ERROR:test: with_suffix: " // p2%path()
call p1%copy_file(p2%path(), overwrite=.true.)
if (.not. p2%is_file()) error stop "ERROR:test: copy_file"

!> same_file, remove
p2 = path_t("./touch.txt")
if (.not. p1%same_file(p2)) error stop "ERROR:test: same_file"
call p2%remove()

!> is_reserved, is_char_device
p1 = path_t(devnull())
b = p1%is_reserved()
b = p1%is_char_device()

!> mkdir
p1 = path_t("test_dir")
call p1%mkdir(ok=b)
if (.not. b) error stop "ERROR:test: mkdir"

!> is_symlink, create_symlink, read_symlink
p1 = path_t("touch.txt")
call p1%touch()

p2 = path_t("touch.link")
if(p2%is_file()) call p2%remove()

call p1%create_symlink(p2%path())
if (.not. p2%is_symlink()) error stop "ERROR:test: create_symlink: " // p2%path()

if (.not. is_windows()) then
p2 = p2%read_symlink()
if (p2%path() /= p1%path()) error stop "ERROR:test: read_symlink " // p2%path()
end if

!> cleanup
call p2%remove()


end block valgrind

print '(a)', "OK: path methods"

end program

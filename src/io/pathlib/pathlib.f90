module pathlib

use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

implicit none (type, external)
private
public :: mkdir, copyfile, expanduser, home, suffix, &
filesep_windows, filesep_unix, &
directory_exists, assert_directory_exists, assert_file_exists, &
make_absolute, is_absolute, get_filename, parent, file_name, stem

interface  ! pathlib_{unix,windows}.f90
module subroutine copyfile(source, dest)
character(*), intent(in) :: source, dest
end subroutine copyfile

module subroutine mkdir(path)
character(*), intent(in) :: path
end subroutine mkdir

module logical function is_absolute(path)
character(*), intent(in) :: path
end function is_absolute

end interface

interface !< pathlib_{intel,gcc}.f90
module logical function directory_exists(path) result(exists)
character(*), intent(in) :: path
end function directory_exists
end interface


contains

pure function suffix(filename)
!! extracts path suffix, including the final "." dot
character(*), intent(in) :: filename
character(:), allocatable :: suffix

integer :: i

i = index(filename, '.', back=.true.)

if (i > 1) then
  suffix = trim(filename(i:))
else
  suffix = ''
end if

end function suffix


pure function parent(instr)

character(*), intent(in) :: instr
character(:), allocatable :: parent

character(len(instr)) :: work
integer :: i

work = filesep_unix(instr)

i = index(work, "/", back=.true.)
if (i > 0) then
  parent = work(:i-1)
else
  parent = "."
end if

end function parent


pure function file_name(instr)

character(*), intent(in) :: instr
character(:), allocatable :: file_name

character(len(instr)) :: work

work = filesep_unix(instr)

file_name = trim(work(index(work, "/", back=.true.) + 1:))

end function file_name


pure function stem(instr)

character(*), intent(in) :: instr
character(:), allocatable :: stem

character(len(instr)) :: work
integer :: i

work = file_name(instr)

i = index(work, '.', back=.true.)
if (i > 0) then
  stem = work(:i - 1)
else
  stem = work
endif

end function stem


function get_filename(path, fstem) result(fn)
!! given a path and stem, find the full filename
!! assumes:
!! 1. "stem" is the file name we wish to find (without suffix or directories)
!! 2. a file exists with suffix (else error)
character(*), intent(in) :: path
character(*), intent(in), optional :: fstem

character(:), allocatable :: fn, path1
integer :: i, L
logical :: exists
character(*), parameter :: suff(3) = [character(4) :: '.h5', '.nc', '.dat']

fn = trim(path)  !< first to avoid undefined return

if(len(fn) == 0) return

if(present(fstem)) then
  if(index(fn, fstem, back=.true.) == 0) then
    !> assume we wish to append stem to path
    fn = fn // '/' // fstem
  elseif(index(fn, '.', back=.true.) > 4) then
    !> it's a stem-matching full path with a suffix
    inquire(file=fn, exist=exists)
    if(.not. exists) fn = ''
    return
  endif
endif

inquire(file=fn, exist=exists)
if(exists) return

path1 = fn

do i = 1, size(suff)
  fn = path1 // trim(suff(i))
  inquire(file=fn, exist=exists)
  if (exists) return
enddo

fn = ''
if(present(fstem)) then
  write(stderr,*) 'ERROR:pathlib:get_filename: ',fstem,' not found in ', path
else
  write(stderr,*) 'ERROR:pathlib:get_filename: file not found: ',path
endif

end function get_filename


function make_absolute(path, top_path) result(abspath)
!! if path is absolute, return expanded path
!! if path is relative, top_path / path

!! NOTE: can only allocate once when it's a function, it will ignore later allocates
!! NOTE: only idempotent if top_path is absolute

character(:), allocatable :: abspath, p
logical :: exists, is_abs
character(*), intent(in) :: path, top_path

p = expanduser(path)
if (is_absolute(p)) then
  abspath = p
else
  abspath = expanduser(top_path) // '/' // p
endif

end function make_absolute


subroutine assert_directory_exists(path)
!! throw error if directory does not exist
character(*), intent(in) :: path

if (.not. directory_exists(path)) error stop 'directory does not exist ' // path

end subroutine assert_directory_exists


subroutine assert_file_exists(path)
!! throw error if file does not exist

character(*), intent(in) :: path
logical :: exists

inquire(file=expanduser(path), exist=exists)

if (.not. exists) error stop 'file does not exist ' // path

end subroutine assert_file_exists


pure function filesep_windows(path) result(swapped)
!! '/' => '\' for Windows systems

character(*), intent(in) :: path
character(len(path)) :: swapped
integer :: i

swapped = path
i = index(swapped, '/')
do while (i > 0)
  swapped(i:i) = char(92)
  i = index(swapped, '/')
end do

end function filesep_windows


pure function filesep_unix(path) result(swapped)
!! '\' => '/'

character(*), intent(in) :: path
character(len(path)) :: swapped
integer :: i

swapped = path
i = index(swapped, char(92))
do while (i > 0)
  swapped(i:i) = '/'
  i = index(swapped, char(92))
end do

end function filesep_unix


function expanduser(in) result (out)
!! resolve home directory as Fortran does not understand tilde
!! works for Linux, Mac, Windows, etc.
character(:), allocatable :: out, homedir
character(*), intent(in) :: in

out = filesep_unix(in)

if (len_trim(out) < 1 .or. out(1:1) /= '~') then
  !! nothing to expand
  out = trim(adjustl(out))
  return
endif

homedir = home()
if (len_trim(homedir) == 0) then
  !! could not determine the home directory
  out = trim(adjustl(out))
  return
endif

if (len_trim(out) < 3) then
  !! ~ or ~/
  out = homedir
else
  !! ~/...
  out = homedir // trim(adjustl(out(3:)))
endif

end function expanduser


function home()
!! https://en.wikipedia.org/wiki/Home_directory#Default_home_directory_per_operating_system
character(:), allocatable :: home
character(256) :: buf
integer :: L, istat

call get_environment_variable("HOME", buf, length=L, status=istat)
if (L==0 .or. istat /= 0) then
  call get_environment_variable("USERPROFILE", buf, length=L, status=istat)
endif

if (L==0 .or. istat /= 0) then
  write(stderr,*) 'ERROR: could not determine home directory from env variable'
  if (istat==1) write(stderr,*) 'env variable does not exist.'
  home = ""
else
  home = trim(buf) // '/'
endif

end function home

end module pathlib

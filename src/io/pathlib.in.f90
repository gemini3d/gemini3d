module pathlib

use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

implicit none (type, external)
private
public :: mkdir, copyfile, expanduser, home, get_suffix, filesep_swap, &
  assert_directory_exists, assert_file_exists,&
  make_absolute, is_absolute, get_filename, parent, file_name

interface  ! pathlib_{unix,windows}.f90
module integer function copyfile(source, dest) result(istat)
character(*), intent(in) :: source, dest
end function copyfile

module integer function mkdir(path) result(istat)
character(*), intent(in) :: path
end function mkdir

module logical function is_absolute(path)
character(*), intent(in) :: path
end function is_absolute

end interface


contains

pure function get_suffix(filename)
!! extracts path suffix, including the final "." dot
character(*), intent(in) :: filename
character(:), allocatable :: get_suffix

get_suffix = filename(index(filename, '.', back=.true.) : len(filename))

end function get_suffix


pure function parent(instr, filesep)

character(*), intent(in) :: instr
character(1), intent(in), optional :: filesep
character(:), allocatable :: parent

character(1) :: sep
integer :: i

sep = '/'
if(present(filesep)) sep = filesep

i = scan(instr, sep, back=.true.)
parent = instr(1:i-1)

end function parent


pure function file_name(instr, filesep)

character(*), intent(in) :: instr
character(1), intent(in), optional :: filesep
character(:), allocatable :: file_name

character(1) :: sep
integer :: i

sep = '/'
if(present(filesep)) sep = filesep

i = scan(instr, sep, back=.true.)
file_name = instr(i+1:len(instr))

end function file_name


function get_filename(path, stem) result(fn)
!! given a path and stem, find the full filename
!! assumes:
!! 1. "stem" is the file name we wish to find (without suffix or directories)
!! 2. a file exists with suffix (else error)
character(*), intent(in) :: path
character(*), intent(in), optional :: stem

character(:), allocatable :: fn, path1
integer :: i, L
logical :: exists
character(*), parameter :: suffix(3) = [character(4) :: '.h5', '.nc', '.dat']

fn = trim(path)  !< first to avoid undefined return

if(len(fn) == 0) return

if(present(stem)) then
  if(index(fn, stem, back=.true.) == 0) then
    !> assume we wish to append stem to path
    fn = fn // '/' // stem
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

do i = 1, size(suffix)
  fn = path1 // trim(suffix(i))
  inquire(file=fn, exist=exists)
  if (exists) return
enddo

fn = ''
if(present(stem)) then
  write(stderr,*) 'ERROR:pathlib:get_filename: ',stem,' not found in ', path
else
  write(stderr,*) 'ERROR:pathlib:get_filename: file not found: ',path
endif

end function get_filename


function make_absolute(path, top_path) result(abspath)
!! if path is absolute, return expanded path
!! if path is relative, make absolute path under absolute top_path

!! NOTE:
!! 1. can only allocate once when it's a function, it will ignore later allocates
!! 2. need trim(adjustl()) to sanitize fixed length namelist input

character(:), allocatable :: abspath, p, t
logical :: exists, is_abs
character(*), intent(in) :: path, top_path

p = expanduser(path)
if (is_absolute(p)) then
  abspath = p
  return
endif

t = expanduser(top_path)
if (.not. is_absolute(t)) write(stderr,*) "WARNING: make_absolute: top_path is not absolute: " // t
abspath = t // '/' // p

end function make_absolute


subroutine assert_directory_exists(path)
!! throw error if directory does not exist
character(*), intent(in) :: path
logical :: exists

@dir_exist@

if (exists) return

error stop 'directory does not exist ' // path

end subroutine assert_directory_exists


subroutine assert_file_exists(path)
!! throw error if file does not exist

character(*), intent(in) :: path
logical :: exists

inquire(file=path, exist=exists)

if (exists) return

error stop 'ERROR: file does not exist ' // path

end subroutine assert_file_exists


function filesep_swap(path) result(swapped)
!! swaps '/' to '\' for Windows systems

character(*), intent(in) :: path
character(len(path)) :: swapped
integer :: i

swapped = path
do
  i = index(swapped, '/')
  if (i == 0) exit
  swapped(i:i) = char(92)
end do

end function filesep_swap


function expanduser(indir)
!! resolve home directory as Fortran does not understand tilde
!! works for Linux, Mac, Windows, etc.
character(:), allocatable :: expanduser, homedir
character(*), intent(in) :: indir

if (len_trim(indir) < 1 .or. indir(1:1) /= '~') then
  !! nothing to expand
  expanduser = trim(adjustl(indir))
  return
endif

homedir = home()
if (len_trim(homedir) == 0) then
  !! could not determine the home directory
  expanduser = trim(adjustl(indir))
  return
endif

if (len_trim(indir) < 3) then
  !! ~ or ~/
  expanduser = homedir
else
  !! ~/...
  expanduser = homedir // trim(adjustl(indir(3:)))
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

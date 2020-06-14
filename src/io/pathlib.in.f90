module pathlib

use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

implicit none (type, external)
private
public :: mkdir, copyfile, expanduser, home, get_suffix, filesep_swap, &
  assert_directory_exists, directory_exists, assert_file_exists,&
  make_absolute

interface  ! pathlib_{unix,windows}.f90
module integer function copyfile(source, dest) result(istat)
character(*), intent(in) :: source, dest
end function copyfile

module integer function mkdir(path) result(istat)
character(*), intent(in) :: path
end function mkdir
end interface

contains

pure function get_suffix(filename)
!! extracts path suffix, including the final "." dot
character(*), intent(in) :: filename
character(:), allocatable :: get_suffix

get_suffix = filename(index(filename, '.', back=.true.) : len(filename))

end function get_suffix


function make_absolute(path, top_path) result(abspath)
!! if path is absolute, return expanded path
!! if path is relative, make absolute path under absolute top_path
!!
!! NOTE: To constrain the dimensionality of this problem across operating
!! systems, we require that:
!! 1. top_path is known to be absolute
!! 2. top_path exists
!! 3. if relative, path already exists

character(:), allocatable :: abspath, rel, top
logical :: exists
character(*), intent(in) :: path, top_path


rel = expanduser(path)
top = expanduser(top_path)

if (.not.directory_exists(top)) then
  write(stderr,*) 'ERROR: pathlib:make_absolute: ',top,' does not exist'
  error stop
endif

inquire(file=rel, exist=exists)
if (exists) then
  abspath = rel
else
  abspath = top // '/' // rel
endif


end function make_absolute


subroutine assert_directory_exists(path)
!! throw error if directory does not exist
character(*), intent(in) :: path

if (directory_exists(path)) return

write(stderr,*) path // ' directory does not exist'
error stop

end subroutine assert_directory_exists


logical function directory_exists(path) result(exists)
!! check if directory exists, handling Intel compiler-specific behavior
character(*), intent(in) :: path

print *, 'check exist ', path
@dir_exist@

end function directory_exists


subroutine assert_file_exists(path)
  !! throw error if file does not exist
  !! this accomodates non-Fortran 2018 error stop with variable character

character(*), intent(in) :: path
logical :: exists

inquire(file=path, exist=exists)

if (exists) return

write(stderr,'(A)') 'ERROR: file does not exist ' // path
error stop

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

subroutine destructor(self)
type(path_t), intent(inout) :: self
if(allocated(self%path_str)) deallocate(self%path_str)
end subroutine destructor


type(path_t) function set_path(path)
character(*), intent(in) :: path

allocate(character(max_path()) :: set_path%path_str)

set_path%path_str = path
end function set_path


pure function get_path(self, istart, iend)
character(:), allocatable :: get_path
class(path_t), intent(in) :: self
integer, intent(in), optional :: istart, iend
integer :: i1, i2

i1 = 1
i2 = len_trim(self%path_str)
if(present(istart)) i1 = istart
if(present(iend))   i2 = iend

get_path = self%path_str(i1:i2)
end function get_path


pure integer function length(self)
!! returns string length len_trim(path)
class(path_t), intent(in) :: self
length = len_trim(self%path_str)
end function


subroutine m_unlink(self, ok)
!! delete the file
class(path_t), intent(in) :: self
logical, intent(out), optional :: ok
call remove(self%path_str, ok)
end subroutine


subroutine m_rename(self, to, ok)
!! rename file "to"
class(path_t), intent(in) :: self
character(*), intent(in) :: to
logical, intent(out), optional :: ok
call fs_rename(self%path_str, to, ok)
end subroutine


function m_as_posix(self) result(r)
!! force Posix "/" file separator
class(path_t), intent(in) :: self
type(path_t) :: r
r%path_str = as_posix(self%path_str)
end function


function m_canonical(self, strict) result(r)
class(path_t), intent(in) :: self
logical, intent(in), optional :: strict
type(path_t) :: r
r%path_str = canonical(self%path_str, strict)
end function


function m_resolve(self, strict) result(r)
class(path_t), intent(in) :: self
logical, intent(in), optional :: strict
type(path_t) :: r
r%path_str = resolve(self%path_str, strict)
end function


subroutine m_set_permissions(self, readable, writable, executable, ok)
class(path_t), intent(in) :: self
logical, intent(in), optional :: readable, writable, executable
logical, intent(out), optional :: ok
call set_permissions(self%path_str, readable, writable, executable, ok)
end subroutine


character(9) function m_get_permissions(self) result (r)
!! get file permissions as POSIX string
class(path_t), intent(in) :: self
r = get_permissions(self%path_str)
end function


subroutine m_copy_file(self, dest, overwrite, ok)
!! copy file from source to destination
!! OVERWRITES existing destination files
class(path_t), intent(in) :: self
character(*), intent(in) :: dest
logical, intent(in), optional :: overwrite
logical, intent(out), optional :: ok
call copy_file(self%path_str, dest, overwrite, ok)
end subroutine


subroutine m_mkdir(self, ok)
!! create a directory, with parents if needed
class(path_t), intent(in) :: self
logical, intent(out), optional :: ok
call mkdir(self%path_str, ok=ok)
end subroutine


function m_read_symlink(self) result(r)
class(path_t), intent(in) :: self
type(path_t) :: r
r%path_str = read_symlink(self%path_str)
end function


subroutine m_create_symlink(self, link, ok)
class(path_t), intent(in) :: self
character(*), intent(in) :: link
logical, intent(out), optional :: ok
call create_symlink(self%path_str, link, ok)
end subroutine


logical function m_exists(self) result(r)
class(path_t), intent(in) :: self
r = exists(self%path_str)
end function


function m_expanduser(self) result(r)
!! resolve home directory as Fortran does not understand tilde
!! works for Linux, Mac, Windows, etc.
class(path_t), intent(in) :: self
type(path_t) :: r
r%path_str = expanduser(self%path_str)
end function


function m_file_name(self) result (r)
!! returns file name without path
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = file_name(self%path_str)
end function


function m_get_owner_name(self) result (r)
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = get_owner_name(self%path_str)
end function


function m_get_owner_group(self) result(r)
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = get_owner_group(self%path_str)
end function


logical function m_is_safe_name(self) result(r)
!! is path a safe name for this filesystem
class(path_t), intent(in) :: self
r = is_safe_name(self%path_str)
end function


integer(int64) function m_file_size(self) result(r)
class(path_t), intent(in) :: self
r = file_size(self%path_str)
end function


logical function m_is_absolute(self) result(r)
!! is path absolute
!! do NOT expanduser() to be consistent with Python etc. filesystem
class(path_t), intent(in) :: self
r = is_absolute(self%path_str)
end function


logical function m_is_char_device(self) result(r)
class(path_t), intent(in) :: self
r = is_char_device(self%path_str)
end function


logical function m_is_dir(self) result(r)
class(path_t), intent(in) :: self
r = is_dir(self%path_str)
end function


logical function m_is_exe(self) result(r)
class(path_t), intent(in) :: self
r = is_exe(self%path_str)
end function


logical function m_is_file(self) result(r)
class(path_t), intent(in) :: self
r = is_file(self%path_str)
end function


logical function m_is_readable(self) result(r)
class(path_t), intent(in) :: self
r = is_readable(self%path_str)
end function


logical function m_is_writable(self) result(r)
class(path_t), intent(in) :: self
r = is_writable(self%path_str)
end function


logical function m_is_reserved(self) result(r)
class(path_t), intent(in) :: self
r = is_reserved(self%path_str)
end function


logical function m_is_symlink(self) result(r)
class(path_t), intent(in) :: self
r = is_symlink(self%path_str)
end function


function m_join(self, other) result(r)
!! returns path_t object with other appended to self using posix separator
type(path_t) :: r
class(path_t), intent(in) :: self
character(*), intent(in) :: other
r%path_str = join(self%path_str, other)
end function


logical function m_is_prefix(self, path) result(r)
class(path_t), intent(in) :: self
character(*), intent(in) :: path
r = is_prefix(self%path_str, path)
end function



logical function m_is_subdir(self, dir) result(r)
class(path_t), intent(in) :: self
character(*), intent(in) :: dir
r = is_subdir(self%path_str, dir)
end function


function m_parent(self) result(r)
!! returns parent directory of path
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = parent(self%path_str)
end function


function m_normal(self) result(r)
!! lexically normalize path
class(path_t), intent(in) :: self
type(path_t) :: r
r%path_str = normal(self%path_str)
end function


function m_relative_to(self, other) result(r)
!! returns other relative to self
class(path_t), intent(in) :: self
character(*), intent(in) :: other
character(:), allocatable :: r

r = relative_to(self%path_str, other)
end function


function m_proximate_to(self, other) result(r)
!! returns other proximate to self
class(path_t), intent(in) :: self
character(*), intent(in) :: other
character(:), allocatable :: r

r = proximate_to(self%path_str, other)
end function


function m_root(self) result(r)
!! returns root of path
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = root(self%path_str)
end function

function m_root_name(self) result(r)
!! returns root name of path
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = root_name(self%path_str)
end function


logical function m_same_file(self, other) result(r)
class(path_t), intent(in) :: self, other
r = same_file(self%path_str, other%path_str)
end function


function m_stem(self) result(r)
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = stem(self%path_str)
end function


function m_suffix(self) result(r)
!! extracts path suffix, including the final "." dot
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = suffix(self%path_str)
end function


subroutine m_touch(self, ok)
class(path_t), intent(in) :: self
logical, intent(out), optional :: ok
call touch(self%path_str, ok)
end subroutine


integer(C_LONG_LONG) function m_get_modtime(self) result(r)
class(path_t), intent(in) :: self
r = get_modtime(self%path_str)
end function


logical function m_set_modtime(self) result(r)
class(path_t), intent(in) :: self
r = set_modtime(self%path_str)
end function


function m_with_suffix(self, new) result(r)
!! replace file suffix with new suffix
class(path_t), intent(in) :: self
type(path_t) :: r
character(*), intent(in) :: new
r%path_str = with_suffix(self%path_str, new)
end function


function m_shortname(self) result(r)
class(path_t), intent(in) :: self
type(path_t) :: r
r%path_str = shortname(self%path_str)
end function


function m_longname(self) result(r)
class(path_t), intent(in) :: self
type(path_t) :: r
r%path_str = longname(self%path_str)
end function


function m_filesystem_type(self) result(r)
class(path_t), intent(in) :: self
character(:), allocatable :: r
r = filesystem_type(self%path_str)
end function

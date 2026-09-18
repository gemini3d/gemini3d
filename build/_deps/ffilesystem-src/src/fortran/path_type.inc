
public :: path_t  !< base class

type :: path_t

private

character(:), allocatable :: path_str

contains

procedure, public :: path=>get_path
procedure, public :: length
procedure, public :: as_posix=>m_as_posix
procedure, public :: join=>m_join
procedure, public :: relative_to=>m_relative_to
procedure, public :: proximate_to=>m_proximate_to
procedure, public :: normal=>m_normal
procedure, public :: exists=>m_exists
procedure, public :: is_char_device=>m_is_char_device
procedure, public :: is_file=>m_is_file
procedure, public :: is_exe=>m_is_exe
procedure, public :: is_readable=>m_is_readable
procedure, public :: is_writable=>m_is_writable
procedure, public :: is_dir=>m_is_dir
procedure, public :: is_safe_name=>m_is_safe_name
procedure, public :: is_prefix=>m_is_prefix
procedure, public :: is_subdir=>m_is_subdir
procedure, public :: is_reserved=>m_is_reserved
procedure, public :: is_absolute=>m_is_absolute
procedure, public :: is_symlink=>m_is_symlink
procedure, public :: read_symlink=>m_read_symlink
procedure, public :: create_symlink=>m_create_symlink
procedure, public :: copy_file=>m_copy_file
procedure, public :: mkdir=>m_mkdir
procedure, public :: touch=>m_touch
procedure, public :: modtime=>m_get_modtime
procedure, public :: set_modtime=>m_set_modtime
procedure, public :: parent=>m_parent
procedure, public :: file_name=>m_file_name
procedure, public :: stem=>m_stem
procedure, public :: root=>m_root
procedure, public :: root_name=>m_root_name
procedure, public :: suffix=>m_suffix
procedure, public :: expanduser=>m_expanduser
procedure, public :: with_suffix=>m_with_suffix
procedure, public :: canonical=>m_canonical
procedure, public :: resolve=>m_resolve
procedure, public :: same_file=>m_same_file
procedure, public :: remove=>m_unlink
procedure, public :: rename=>m_rename
procedure, public :: file_size=>m_file_size
procedure, public :: set_permissions=>m_set_permissions
procedure, public :: get_permissions=>m_get_permissions
procedure, public :: shortname=>m_shortname
procedure, public :: longname=>m_longname
procedure, public :: filesystem_type=>m_filesystem_type
procedure, public :: get_owner_name=>m_get_owner_name
procedure, public :: get_owner_group=>m_get_owner_group

final :: destructor

end type path_t

interface path_t
!! constructor
  module procedure set_path
end interface

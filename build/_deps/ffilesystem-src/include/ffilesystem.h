#ifndef FFILESYSTEM_H
#define FFILESYSTEM_H

#include <sys/types.h> // for pid_t, ino_t, dev_t, mode_t, ...
// https://www.man7.org/linux/man-pages/man3/stat.3type.html

#ifdef _MSC_VER

#ifndef pid_t
#ifdef __cplusplus
using pid_t = long;
#else
typedef long pid_t;
#endif
#endif

#ifndef mode_t
#ifdef __cplusplus
using mode_t = unsigned int;
#else
typedef unsigned int mode_t;
#endif
#endif

#endif // _MSC_VER

#if !defined(__APPLE__) && (defined(BSD) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__))
#define FFS_BSD 1
#endif
#if defined(__APPLE__) && defined(__MACH__)
#define FFS_DARWIN 1
#endif


#ifdef __cplusplus

#include <cstdint> // for uintmax_t

constexpr int fs_trace =
#ifdef FS_TRACE
  1;
#else
  0;
#endif

constexpr std::uintmax_t fs_unknown_size = static_cast<std::uintmax_t>(-1);

#if !defined(__has_cpp_attribute)
// this is a C++20 feature, but available in many older compilers
#  define __has_cpp_attribute(x)  0
#endif

// GCC itself does it this way https://github.com/gcc-mirror/gcc/blob/78b56a12dd028b9b4051422c6bad6260055e4465/libcpp/system.h#L426
#if __has_cpp_attribute(unlikely)
#  define FFS_UNLIKELY [[unlikely]]
#  define FFS_LIKELY [[likely]]
#else
#  define FFS_UNLIKELY
#  define FFS_LIKELY
#endif

#include <cstdlib>
#include <vector>
#include <string>
#include <string_view>
#include <ctime> // for time_t
#include <system_error>
#include <optional>

#if __has_include(<source_location>)
#include <source_location>
#endif

#ifdef HAVE_CXX_FILESYSTEM

#include <filesystem>

std::optional<std::filesystem::file_time_type> fs_get_modtime_fs(std::string_view);

#endif

std::string fs_absolute(std::string_view);
std::string fs_absolute(std::string_view, std::string_view);

void fs_as_posix(std::string&);
std::string fs_as_posix(std::string_view);

void fs_as_windows(std::string&);

std::string fs_backend();

std::string fs_canonical(std::string_view, const bool = false);

std::string fs_compiler();

bool fs_copy_file(std::string_view, std::string_view, bool);

std::string fs_cpu_arch();

bool fs_create_symlink(std::string_view, std::string_view);

std::string fs_drop_slash(std::string_view);
void fs_drop_trailing_slash(std::string&);

bool fs_equivalent(std::string_view, std::string_view);

bool fs_is_removable(std::string_view);
bool fs_exists(std::string_view);
bool fs_lexists(std::string_view);

std::string fs_exe_path();

std::string fs_expanduser(std::string_view);

std::string fs_file_name(std::string_view);

std::uintmax_t fs_file_size(std::string_view);

std::string fs_filesystem_type(std::string_view);

std::string fs_generate_random_alphanumeric_string(const std::string::size_type);
void fs_ascii_lower(std::string&);
void fs_ascii_upper(std::string&);

std::string fs_get_cwd();

std::string fs_get_homedir();

std::string fs_get_locale_name();

std::time_t fs_get_modtime(std::string_view);

std::string fs_get_owner_group(std::string_view);
std::string fs_get_owner_name(std::string_view);

std::string fs_get_permissions(std::string_view);

std::string fs_get_profile_dir();

std::string fs_get_shell();

std::string fs_get_tempdir();

std::string fs_get_terminal();

std::string fs_get_username();

std::optional<std::string> fs_getenv(std::string_view);

std::uintmax_t fs_hard_link_count(std::string_view);

std::string fs_hostname();
bool fs_has_filename(std::string_view);

bool fs_is_absolute(std::string_view);
bool fs_is_appexec_alias(std::string_view);
bool fs_is_case_sensitive(std::string_view);
bool fs_is_char_device(std::string_view);
bool fs_is_dir(std::string_view);
bool fs_is_empty(std::string_view);

bool fs_is_exe(std::string_view);
bool fs_is_executable_binary(std::string_view);

bool fs_is_fifo(std::string_view);
bool fs_is_file(std::string_view);
bool fs_is_other(std::string_view);

bool fs_is_prefix(std::string_view, std::string_view);
bool fs_is_readable(std::string_view);
bool fs_is_reserved(std::string_view);
bool fs_is_safe_name(std::string_view);
bool fs_is_subdir(std::string_view, std::string_view);
bool fs_is_symlink(std::string_view);
bool fs_is_writable(std::string_view);

std::string fs_join(std::string_view, std::string_view);

std::string fs_lexically_normal(std::string_view);

std::string fs_lib_path();

std::string fs_libc();
std::string fs_libcxx();

std::string fs_longname(std::string_view);

std::string fs_make_preferred(std::string_view);

std::string::size_type fs_max_component(std::string_view);

bool fs_mkdir(std::string_view);

bool fs_non_ascii(std::string_view);

std::string fs_normal(std::string_view);
std::vector<std::string> fs_normal_vector(std::string_view);

std::string fs_os_version();
std::string fs_parent(std::string_view);

void fs_emit_error();
void fs_print_error(std::string_view
#if defined(__cpp_lib_source_location)
, const std::source_location& = std::source_location::current()
#endif
);

void fs_print_error(std::string_view, const std::error_code&
#if defined(__cpp_lib_source_location)
, const std::source_location& = std::source_location::current()
#endif
);

void fs_print_error(std::string_view, std::string_view
#if defined(__cpp_lib_source_location)
, const std::source_location& = std::source_location::current()
#endif
);

void fs_print_error(std::string_view, std::string_view, const std::error_code&
#if defined(__cpp_lib_source_location)
, const std::source_location& = std::source_location::current()
#endif
);

std::string fs_proximate_to(std::string_view, std::string_view);

std::string fs_read_symlink(std::string_view);

std::string fs_realpath(std::string_view);

std::string fs_relative_to(std::string_view, std::string_view);

bool fs_remove(std::string_view);

bool fs_rename(std::string_view, std::string_view);

std::string fs_resolve(std::string_view, const bool = false);

std::string fs_root(std::string_view);
std::string fs_root_name([[maybe_unused]] std::string_view);

bool fs_set_cwd(std::string_view);

bool fs_set_modtime(std::string_view, const bool = false);

bool fs_set_permissions(std::string_view, int, int, int);

bool fs_setenv(std::string_view, std::string_view);

std::string fs_shortname(std::string_view);

bool fs_slash_first(std::string_view);

std::uintmax_t fs_space_available(std::string_view);
std::uintmax_t fs_space_capacity(std::string_view);
std::size_t fs_get_blksize(std::string_view);

std::vector<std::string> fs_split(std::string_view);
std::vector<std::string> fs_split_pathsep(std::string_view);
bool fs_stdin_tty();

std::string fs_stem(std::string_view);

std::string::size_type fs_str2char(std::string_view, char*, const std::string::size_type);

std::string fs_suffix(std::string_view);

bool fs_touch(std::string_view);

std::string fs_to_cygpath(std::string_view);

std::string fs_to_winpath(std::string_view);

void fs_trim(std::string&);

std::string fs_which(std::string_view, std::string_view = {}, const bool = false);

std::string fs_win32_final_path(std::string_view);
std::string fs_win32_full_name(std::string_view);
std::string fs_win32_to_narrow([[maybe_unused]] std::wstring_view);
std::wstring fs_win32_to_wide(std::string_view);
bool fs_win32_is_symlink(std::string_view);
bool fs_win32_is_ext_path(std::string_view);

std::string fs_with_suffix(std::string_view, std::string_view);

std::string_view::size_type fs_symlink_length([[maybe_unused]] std::string_view);

mode_t fs_st_mode(std::string_view);

dev_t fs_st_dev(std::string_view);
ino_t fs_inode(std::string_view);
struct passwd* fs_getpwuid();

// ---------------------------------------------------------------------------

extern "C" {

#else  // C only
#if __STDC_VERSION__ < 202311L
#include <stdbool.h>
#endif

#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>

// internal use function only for C
size_t fs_strncpy(const char*, char*, const size_t);

#endif  // end C-only


#if !defined(__has_c_attribute)
#  define __has_c_attribute(x)  0
#endif

long fs_cpp_lang();
long fs_cpp_format();
long fs_cpp_ranges();

long fs_c_lang();

long libcxx_release();

size_t fs_get_max_path();
long fs_get_max_open_files();
size_t fs_max_component(const char*);

bool fs_is_optimized();

pid_t fs_getpid();

unsigned long long fs_total_sys_memory();
unsigned long long fs_get_free_memory();

char fs_filesep();
char fs_pathsep();
const char* fs_devnull();

bool fs_is_admin();
bool fs_is_android();
bool fs_is_bsd();
bool fs_is_macos();
bool fs_is_linux();
bool fs_is_unix();
bool fs_is_windows();
int fs_is_wsl();
bool fs_is_mingw();
bool fs_is_clangcl();
bool fs_is_msvc();
bool fs_is_appleclang();

bool fs_win32_long_paths_enabled();
bool fs_is_cygwin();
bool fs_is_rosetta();

bool fs_has_statx();

bool fs_is_empty(const char*);

bool fs_is_safe_name(const char*);
bool fs_is_appexec_alias(const char*);

size_t fs_as_posix(char*, const size_t);
size_t fs_as_windows(char*, const size_t);

size_t fs_normal(const char*, char*, const size_t);

size_t fs_realpath(const char*, char*, const size_t);

size_t fs_join(const char*, const char*, char*, const size_t);
size_t fs_absolute(const char*, const char*, char*, const size_t);

size_t fs_file_name(const char*, char*, const size_t);
size_t fs_stem(const char*, char*, const size_t);
size_t fs_parent(const char*, char*, const size_t);
size_t fs_suffix(const char*, char*, const size_t);

size_t fs_root(const char*, char*, const size_t);
size_t fs_root_name(const char*, char*, const size_t);

size_t fs_with_suffix(const char*, const char*, char*, const size_t);

size_t fs_which(const char*, const char*, const bool, char*, const size_t);

bool fs_is_symlink(const char*);
bool fs_create_symlink(const char*, const char*);
size_t fs_read_symlink(const char*, char*, const size_t);

bool fs_mkdir(const char*);

bool fs_exists(const char*);
bool fs_lexists(const char*);
bool fs_is_removable(const char*);

bool fs_has_filename(const char* path);
bool fs_is_absolute(const char*);
bool fs_is_case_sensitive(const char*);

bool fs_is_char_device(const char*);
bool fs_is_fifo(const char*);

bool fs_is_dir(const char*);
bool fs_is_file(const char*);

bool fs_is_other(const char*);

bool fs_is_exe(const char*);
bool fs_is_executable_binary(const char*);

bool fs_is_readable(const char*);
bool fs_is_writable(const char*);
bool fs_is_reserved(const char*);

bool fs_is_prefix(const char*, const char*);
bool fs_is_subdir(const char*, const char*);

bool fs_set_permissions(const char*, int, int, int);

size_t fs_get_permissions(const char*, char*, const size_t);

bool fs_remove(const char*);
bool fs_rename(const char*, const char*);

size_t fs_canonical(const char*, const bool, char*, const size_t);
size_t fs_resolve(const char*, const bool, char*, const size_t);

bool fs_equivalent(const char*, const char*);
bool fs_copy_file(const char*, const char*, bool);

size_t fs_relative_to(const char*, const char*, char*, const size_t);
size_t fs_proximate_to(const char*, const char*, char*, const size_t);

bool fs_touch(const char*);

time_t fs_get_modtime(const char*);
bool fs_set_modtime(const char*);

size_t fs_get_cwd(char*, const size_t);
bool fs_set_cwd(const char*);

size_t fs_get_homedir(char*, const size_t);
size_t fs_get_profile_dir(char*, const size_t);
size_t fs_get_username(char*, const size_t);

size_t fs_get_owner_name(const char*, char*, const size_t);
size_t fs_get_owner_group(const char*, char*, const size_t);

size_t fs_get_tempdir(char*, const size_t);

size_t fs_expanduser(const char*, char*, const size_t);

uintmax_t fs_file_size(const char*);

uintmax_t fs_space_available(const char*);
uintmax_t fs_space_capacity(const char*);
size_t fs_get_blksize(const char*);

ino_t fs_inode(const char*);
dev_t fs_st_dev(const char*);

uintmax_t fs_hard_link_count(const char*);

size_t fs_backend(char*, const size_t);

size_t fs_exe_path(char*, const size_t);
size_t fs_lib_path(char*, const size_t);

size_t fs_compiler(char*, const size_t);
size_t fs_get_shell(char*, const size_t);
size_t fs_get_terminal(char*, const size_t);

size_t fs_hostname(char*, const size_t);

size_t fs_shortname(const char*, char*, const size_t);
size_t fs_longname(const char*, char*, const size_t);

size_t fs_getenv(const char*, char*, const size_t);
bool fs_setenv(const char*, const char*);

size_t fs_filesystem_type(const char*, char*, const size_t);

void fs_print_error(const char*, const char*);

size_t fs_to_cygpath(const char*, char*, const size_t);
size_t fs_to_winpath(const char*, char*, const size_t);

size_t fs_cpu_arch(char*, const size_t);

#ifdef __cplusplus
}
#endif

#endif

// IWYU pragma: no_include <bits/statx-generic.h>
// IWYU pragma: no_include <linux/stat.h>

#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include "ffilesystem.h"

#include <string_view>
#include <system_error>
#include <iostream>  // IWYU pragma: keep
#include <cstdint> // uintmax_t

// include even if <filesystem> is available
#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <io.h> // _access_s
#else
#include <unistd.h>
#endif

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#endif

#include <sys/types.h>  // IWYU pragma: keep
#include <sys/stat.h>   // IWYU pragma: keep

#if __has_include(<fcntl.h>)
#include <fcntl.h>   // AT_* constants for statx
#endif


#if defined(_WIN32)
static DWORD fs_win32_file_type(std::string_view path){

// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-createfilea
  HANDLE h = CreateFileW(fs_win32_to_wide(path).c_str(),
                         0,
                         FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                         nullptr, OPEN_EXISTING, 0, nullptr);

  if(h == INVALID_HANDLE_VALUE){
    DWORD err = GetLastError();
    switch (err) {
      case ERROR_CANT_ACCESS_FILE: case ERROR_FILE_NOT_FOUND: case ERROR_PATH_NOT_FOUND: case ERROR_SUCCESS:
        return FILE_TYPE_UNKNOWN;
      default:
        fs_print_error(path);
        return FILE_TYPE_UNKNOWN;
    }
  }

// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getfiletype
  DWORD t = GetFileType(h);
  CloseHandle(h);
  return t;
}
#endif


bool fs_has_statx()
{
// https://www.man7.org/linux/man-pages/man2/statx.2.html
#if defined(HAVE_STATX)
  return true;
#else
  return false;
#endif
}


mode_t
fs_st_mode(std::string_view path)
{

  const std::string cpath(path);
#if defined(HAVE_STATX)
// Linux Glibc only
// https://www.gnu.org/software/gnulib/manual/html_node/statx.html
// https://www.man7.org/linux/man-pages/man2/statx.2.html

  struct statx x;
  if (::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT, STATX_MODE, &x) == 0) {
    return x.stx_mode;
  } else if (errno != ENOSYS) {
    return 0;
  }
#endif

// https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/stat-functions
  if (struct stat s; ::stat(cpath.c_str(), &s) == 0)
    return s.st_mode;

  return 0;
}


bool
fs_exists(std::string_view path)
{
  // fs_exists() is true even if path is non-readable, as long as UID has permission to see it.
  // this is like Python pathlib.Path.exists()
  // unlike kwSys:SystemTools:FileExists which uses R_OK instead of F_OK like this project.

  bool ok;
#if defined(HAVE_CXX_FILESYSTEM)
  std::error_code ec;
  ok = (Filesystem::exists(path, ec) && !ec) ||
        (fs_is_msvc() && (fs_is_appexec_alias(path) || fs_is_char_device(path)));
  if (ec && ec != std::errc::no_such_file_or_directory)
    fs_print_error(path, ec);
#else

  const std::string cpath(path);

#if defined(_WIN32)
// https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/access-s-waccess-s
  ok = _access_s(cpath.c_str(), 0) == 0 || fs_is_char_device(path);

  // to use the approach below, need more advanced techniques like
  // https://gitlab.kitware.com/cmake/cmake/-/blob/master/Source/kwsys/SystemTools.cxx#L1408

  // WIN32_FILE_ATTRIBUTE_DATA fad;

  // ok = GetFileAttributesExW(fs_win32_to_wide(path).c_str(), GetFileExInfoStandard, &fad);
  // if (!ok){
  //   DWORD err = GetLastError();
  //   if (err != ERROR_FILE_NOT_FOUND && err != ERROR_PATH_NOT_FOUND)
  //     fs_print_error(path);
  // }
#else
  // https://www.man7.org/linux/man-pages/man2/access.2.html
  ok = access(cpath.c_str(), F_OK) == 0;
  if (!ok && errno != ENOENT && errno != ENOTDIR)
    fs_print_error(path);
#endif

#endif

  return ok;
}


bool
fs_is_dir(std::string_view path)
{
  // is path a directory or a symlink to a directory

  bool ok;
#if defined(HAVE_CXX_FILESYSTEM)
// NOTE: Windows top-level drive "C:" needs a trailing slash "C:/"
  std::error_code ec;
  ok = Filesystem::is_directory(path, ec);
  if (ec && ec != std::errc::no_such_file_or_directory)
    fs_print_error(path, ec);
#elif defined(_WIN32)
// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getfileattributesexa
// this also works for Symlinks to directories
  WIN32_FILE_ATTRIBUTE_DATA fad;

  ok = GetFileAttributesExW(fs_win32_to_wide(path).c_str(), GetFileExInfoStandard, &fad) &&
       (fad.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY);
#else
  ok = S_ISDIR(fs_st_mode(path));
#endif

  return ok;
}


bool
fs_is_file(std::string_view path)
{
  // is path a regular file or a symlink to a regular file.
  // not a directory, device, or symlink to a directory.
  // stat() doesn't detect App Execution Aliases
  // AppExec Alias have FILE_ATTRIBUTE_REPARSE_POINT | FILE_ATTRIBUTE_ARCHIVE
  // but need to check the reparse point data for IO_REPARSE_TAG_APPEXECLINK

  bool ok;
#if defined(HAVE_CXX_FILESYSTEM)
  std::error_code ec;
  ok = (Filesystem::is_regular_file(path, ec) && !ec) ||
        (fs_is_msvc() && fs_is_appexec_alias(path));
#elif defined(_WIN32)
  ok = fs_win32_file_type(path) == FILE_TYPE_DISK || fs_is_appexec_alias(path);
#else
  ok = S_ISREG(fs_st_mode(path));
#endif

  return ok;
}


bool
fs_is_fifo(std::string_view path)
{
  // mkfifo() or CreateNamedPipe()
  bool ok;

#if defined(_WIN32)
  ok = fs_win32_file_type(path) == FILE_TYPE_PIPE;
#elif defined(HAVE_CXX_FILESYSTEM)
  std::error_code ec;
  ok = Filesystem::is_fifo(path, ec) && !ec;
#else
  ok = S_ISFIFO(fs_st_mode(path));
#endif

  return ok;
}


bool fs_is_char_device(std::string_view path)
{
// character device like /dev/null or CONIN$

  bool ok;
#if defined(_WIN32)
// currently broken in MSVC STL and MinGW Clang ARM for <filesystem>
  ok = fs_win32_file_type(path) == FILE_TYPE_CHAR;
#elif defined(HAVE_CXX_FILESYSTEM)
  std::error_code ec;
  ok = Filesystem::is_character_file(path, ec) && !ec;
#else
  // Windows: https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/fstat-fstat32-fstat64-fstati64-fstat32i64-fstat64i32
  ok = S_ISCHR(fs_st_mode(path));
#endif

  return ok;
}


bool fs_is_other(std::string_view path)
{
  // just as defined by <filesystem>
  // note that the symlink could point to something not a file or directory
  // https://en.cppreference.com/w/cpp/filesystem/is_other
  return fs_exists(path) && !fs_is_file(path) && !fs_is_dir(path) && !fs_is_symlink(path);
}


bool fs_is_readable(std::string_view path)
{
  // is path readable by the user
  // does not guarantee that the path can be opened (for example, it may be locked)
  const std::string cpath(path);

#if defined(_WIN32)
  // MSVC / MinGW ::perms doesn't detect App Execution Aliases readability
  // like Python os.access(path, os.R_OK) or uv_is_readable().
  // same reasons as our fs_is_writable().
  // https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/access-s-waccess-s

  return _access_s(cpath.c_str(), 4) == 0;
#else
  return access(cpath.c_str(), R_OK) == 0;
#endif
}


bool fs_is_writable(std::string_view path)
{
  // is path writable by the user.
  // this is a more strict test than checking permissions bits,
  // because it also checks ACLs and parent directory writability for creating new files.
  // checks that path is accessible, unlink std::filesystem::perms -- ours is a stricter test
  // more in accord with user plain expectations of "writable"
  // and with Python's os.access(path, os.W_OK) or uv_is_writable().
  // std::filesystem::perms are not as useful because they don't check ACLs
  // or other platform-specific permissions, and they don't check writability of parent directories
  // for creating new files.

  std::string cpath(path);

#if defined(_WIN32)
  // https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/access-s-waccess-s
  fs_as_windows(cpath);
  return _access_s(cpath.c_str(), 2) == 0;
#else
  return access(cpath.c_str(), W_OK) == 0;
#endif
}


std::uintmax_t fs_hard_link_count(std::string_view path)
{
  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)

  auto s = Filesystem::hard_link_count(path, ec);
  if(ec)
    fs_print_error(path, ec);

  return s;

#else

  int r = 0;
  const std::string cpath(path);

#if defined(HAVE_STATX)
// https://www.man7.org/linux/man-pages/man2/statx.2.html
  struct statx sx;
  r = ::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT, STATX_NLINK, &sx);
  if (r == 0)
    return sx.stx_nlink;
#endif

  if (r == 0 || errno == ENOSYS){
    if (struct stat s; !::stat(cpath.c_str(), &s))
      return s.st_nlink;
  }

  fs_print_error(path, ec);
  return fs_unknown_size;
#endif
}

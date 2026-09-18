// IWYU pragma: no_include <bits/statx-generic.h>
// IWYU pragma: no_include <linux/stat.h>
// IWYU pragma: no_include <bits/chrono.h>

#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include "ffilesystem.h"

#include <iostream>
#include <optional>
#include <system_error>

#ifdef HAVE_CLOCK_CAST
// this include can break even when functions not used if C++ ABI problem e.g. Clang 14 with GCC 13
#include <chrono>
#endif
#include <ctime>                // for time_t

// always include to avoid dev confusion
#include <sys/types.h> // IWYU pragma: keep
#include <sys/stat.h> // statx, stat

#include <string_view>

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#elif defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#endif

// standalone #if
#if __has_include(<fcntl.h>)
#include <fcntl.h> // utimensat, AT_* constants
#endif


std::time_t fs_get_modtime(std::string_view path)
{

#if defined(HAVE_CLOCK_CAST) && defined(HAVE_CXX_FILESYSTEM)
  if(const auto &t_fs = fs_get_modtime_fs(path); t_fs){
    const auto t_sys = std::chrono::clock_cast<std::chrono::system_clock>(t_fs.value());
    return std::chrono::system_clock::to_time_t(t_sys);
  }
#else

  int r = 0;
  const std::string cpath(path);

#if defined(HAVE_STATX)
// https://www.man7.org/linux/man-pages/man2/statx.2.html
  struct statx sx;
  r = ::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT, STATX_MTIME, &sx);
  if (r == 0)
    return sx.stx_mtime.tv_sec;
#endif
  if (r == 0 || errno == ENOSYS){
    if (struct stat s; !::stat(cpath.c_str(), &s))
      return s.st_mtime;
  }
#endif

  fs_print_error(path);
  return {};
}


#ifdef HAVE_CXX_FILESYSTEM
std::optional<Filesystem::file_time_type> fs_get_modtime_fs(std::string_view path)
{
  std::error_code ec;

  if(Filesystem::file_time_type t_fs = Filesystem::last_write_time(path, ec); !ec)
    return t_fs;

  fs_print_error(path, ec);
  return {};
}
#endif


bool fs_set_modtime(std::string_view path, const bool quiet)
{

  std::error_code ec;

#ifdef HAVE_CXX_FILESYSTEM

  if(Filesystem::last_write_time(path, Filesystem::file_time_type::clock::now(), ec); !ec)
    return true;
  // techinically IWYU <chrono> but that can break some compilers, and it works without the include.

#elif defined(_WIN32)
// https://learn.microsoft.com/en-us/windows/win32/SysInfo/changing-a-file-time-to-the-current-time
  HANDLE h = CreateFileW(fs_win32_to_wide(path).c_str(),
                         FILE_WRITE_ATTRIBUTES,
                         FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                         nullptr,
                         OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS, nullptr);
  if (h != INVALID_HANDLE_VALUE){
    FILETIME t;
    // https://learn.microsoft.com/en-us/windows/win32/api/sysinfoapi/nf-sysinfoapi-getsystemtimeasfiletime
    GetSystemTimeAsFileTime(&t);
    BOOL ok = SetFileTime(h, nullptr, nullptr, &t);

    CloseHandle(h);
    if(ok)
      return true;
  }
#else
  // utimensat available in macOS >= 10.13
  // https://github.com/python/cpython/issues/75782
  // https://gitlab.kitware.com/cmake/cmake/-/issues/17101
  const std::string cpath(path);
  if (::utimensat(AT_FDCWD, cpath.c_str(), nullptr, 0) == 0)
    return true;
#endif

  if (!quiet)
    fs_print_error(path, ec);

  return false;
}

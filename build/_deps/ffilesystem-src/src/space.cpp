#include "ffilesystem.h"

#include <cstdint> // uintmax_t
#include <system_error> // for error_code

#include <string_view>


#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#else

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <Windows.h> // GetDiskFreeSpaceEx
#else

#include <unistd.h>

#if __has_include(<sys/statvfs.h>)
#define HAVE_STATVFS
#include <sys/statvfs.h>
#endif

#endif

#endif


std::uintmax_t fs_space_available(std::string_view path)
{
  // filesystem space available for device holding path

  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)
  auto s = Filesystem::space(path, ec);
  if (ec)
    fs_print_error(path, ec);

  if (ec || s.available == fs_unknown_size)
    return fs_unknown_size;

  return s.available;
#elif defined(_WIN32)
  // https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getdiskfreespaceexw
  if(ULARGE_INTEGER b; GetDiskFreeSpaceExW(fs_win32_to_wide(path).c_str(), &b, nullptr, nullptr) != 0)
    return b.QuadPart;
#elif defined(HAVE_STATVFS)
  // https://www.man7.org/linux/man-pages/man3/statvfs.3.html
  // https://unix.stackexchange.com/a/703650
  struct statvfs stat;
  const std::string cpath(path);

  if (!::statvfs(cpath.c_str(), &stat))
    return (stat.f_frsize ? stat.f_frsize : stat.f_bsize) * stat.f_bavail;
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(path, ec);
  return fs_unknown_size;
}


std::uintmax_t fs_space_capacity(std::string_view path)
{
  // filesystem space capacity for device holding path
  // total size of the filesystem, in bytes
  // This is the quantity available to the non-privileged user,
  // not the total physical disk size.

  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)
  auto s = Filesystem::space(path, ec);
  if (ec)
    fs_print_error(path, ec);

  if (ec || s.capacity == fs_unknown_size)
    return fs_unknown_size;

  return s.capacity;
#elif defined(_WIN32)
  if(ULARGE_INTEGER b; GetDiskFreeSpaceExW(fs_win32_to_wide(path).c_str(), nullptr, &b, nullptr) != 0)
    return b.QuadPart;
#elif defined(HAVE_STATVFS)
  struct statvfs stat;
  const std::string cpath(path);

  if (!::statvfs(cpath.c_str(), &stat))
    return (stat.f_frsize ? stat.f_frsize : stat.f_bsize) * stat.f_blocks;
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(path, ec);
  return fs_unknown_size;
}

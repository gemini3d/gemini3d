#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include <string>
#include <system_error>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <winioctl.h> // DeviceIoControl
#endif

#include <cerrno>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(HAVE_STATX)
#include <sys/sysmacros.h> // for makedev
#endif

#if __has_include(<fcntl.h>)
#include <fcntl.h>   // AT_* constants for statx
#endif

#include "ffilesystem.h"


std::size_t fs_get_blksize(std::string_view path)
{
  // block size in bytes
#if defined(_WIN32)

  const std::string root = fs_root_name(path);
  if (root.empty())
    return {};

  HANDLE h = CreateFileW(fs_win32_to_wide(R"(\\.\)" + root).c_str(),
                         0,
                         FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                         nullptr, OPEN_EXISTING, 0, nullptr);

  if (h != INVALID_HANDLE_VALUE) {
    DISK_GEOMETRY_EX diskGeometry = {};
    ZeroMemory(&diskGeometry, sizeof(DISK_GEOMETRY_EX));
    DWORD bytesReturned = 0;

    bool const ok = DeviceIoControl(h, IOCTL_DISK_GET_DRIVE_GEOMETRY_EX,
                        nullptr, 0, &diskGeometry, sizeof(diskGeometry),
                        &bytesReturned, nullptr);

    CloseHandle(h);
    if (ok)
      return diskGeometry.Geometry.BytesPerSector;
  }

#else

  const std::string cpath(path);
  int r = 0;

#if defined(HAVE_STATX)
  struct statx sx;
  r = ::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT | AT_SYMLINK_NOFOLLOW, STATX_BASIC_STATS, &sx);
  if (r == 0)
    return sx.stx_blksize;
#endif

  if (r == 0 || errno == ENOSYS){
    if (struct stat s; !::stat(cpath.c_str(), &s))
      return s.st_blksize;
  }
#endif

  fs_print_error(path);
  return {};
}


dev_t fs_st_dev(std::string_view path)
{
  // device number of the file or directory
  // in general dev_t may be unsigned.

  const std::string cpath(path);
  int r = 0;

#if defined(HAVE_STATX)

  struct statx x;
  r = ::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT, STATX_INO, &x);

  if (r == 0)
    return makedev(x.stx_dev_major, x.stx_dev_minor);
  // don't call as ::makedev because some platforms e.g. Android have makedev as a macro

#endif

  if (r == 0 || errno == ENOSYS){
    if(struct stat s; ::stat(cpath.c_str(), &s) == 0)
      return s.st_dev;
  }

  fs_print_error(path);
  return {};
}


ino_t fs_inode(std::string_view path)
{
  // inode number of the file or directory
  //
  // Windows: .st_ino is always zero.
  // See source code for fs_equivalent() for how to use BY_HANDLE_FILE_INFORMATION
  // with GetFileInformationByHandle().

  std::error_code ec;

#if defined(_WIN32)

#if defined(HAVE_GETFILEINFORMATIONBYNAME)
  // https://learn.microsoft.com/en-us/windows-hardware/drivers/ddi/ntifs/ns-ntifs-file_stat_basic_information
  FILE_STAT_BASIC_INFORMATION f1;

  if (GetFileInformationByName(fs_win32_to_wide(path).c_str(), FileStatBasicByNameInfo, &f1, sizeof(f1)))
    return static_cast<ino_t>(f1.FileId.QuadPart);

#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

#else

  const std::string cpath(path);
  int r = 0;

#if defined(HAVE_STATX)

  struct statx x;
  r = ::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT, STATX_INO, &x);
  if (r == 0)
    return x.stx_ino;

#endif

  if (r == 0 || errno == ENOSYS) {
    if(struct stat s; ::stat(cpath.c_str(), &s) == 0)
      return s.st_ino;
  }

#endif // _WIN32

  fs_print_error(path, ec);
  return 0;

}

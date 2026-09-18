#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include "ffilesystem.h"

#include <string_view>
#include <system_error>
#include <iostream>  // IWYU pragma: keep

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#else
#include <sys/types.h>
#include <sys/stat.h>

#if __has_include(<fcntl.h>)
#include <fcntl.h>   // AT_* constants for statx
#endif

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#endif

#endif


#if !defined(HAVE_CXX_FILESYSTEM) && defined(_WIN32)
static bool fs_win32_equiv(std::string_view path1, std::string_view path2)
{
// GetFileInformationByName in Windows >= 24H2
// https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-getfileinformationbyname
// https://github.com/rust-lang/rust/issues/130169

std::error_code ec;

#if defined(HAVE_GETFILEINFORMATIONBYNAME)
  // https://learn.microsoft.com/en-us/windows-hardware/drivers/ddi/ntifs/ns-ntifs-file_stat_basic_information
  FILE_STAT_BASIC_INFORMATION f1, f2;

 if ( GetFileInformationByName(fs_win32_to_wide(path1).c_str(), FileStatBasicByNameInfo, &f1, sizeof(f1)) &&
      GetFileInformationByName(fs_win32_to_wide(path2).c_str(), FileStatBasicByNameInfo, &f2, sizeof(f2))) {
        return f1.VolumeSerialNumber.QuadPart == f2.VolumeSerialNumber.QuadPart &&
               f1.FileId.QuadPart == f2.FileId.QuadPart;
  }
  // .FileID and .VolumeSerialNumber are LARGE_INTEGER

#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(path1, path2, ec);
  return false;
}
#endif


bool fs_equivalent(std::string_view path1, std::string_view path2)
{
  // non-existent paths are not equivalent

  std::error_code ec;

#ifdef HAVE_CXX_FILESYSTEM

  if(bool e = Filesystem::equivalent(path1, path2, ec); !ec)
    return e;

#else

#if defined(_WIN32)

  return fs_win32_equiv(path1, path2);

#else
  int r1 = 0;
  int r2 = 0;
  const std::string p1(path1);
  const std::string p2(path2);

// https://www.man7.org/linux/man-pages/man7/inode.7.html
#if defined(HAVE_STATX)

  struct statx x1, x2;

  r1 = ::statx(AT_FDCWD, p1.c_str(), AT_NO_AUTOMOUNT, STATX_INO, &x1);
  if(r1 == 0){
    r2 = ::statx(AT_FDCWD, p2.c_str(), AT_NO_AUTOMOUNT, STATX_INO, &x2);
    if(r2 == 0)
      return x1.stx_dev_major == x2.stx_dev_major && x1.stx_dev_minor == x2.stx_dev_minor && x1.stx_ino == x2.stx_ino;
  }

#endif

  if((r1 == 0 && r2 == 0) || errno == ENOSYS){
    struct stat s1, s2;

    // https://www.boost.org/doc/libs/1_86_0/libs/filesystem/doc/reference.html#equivalent
    if(!::stat(p1.c_str(), &s1) && !::stat(p2.c_str(), &s2))
      return s1.st_dev == s2.st_dev && s1.st_ino == s2.st_ino;
  }

#endif

#endif // HAVE_CXX_FILESYSTEM

  fs_print_error(path1, path2, ec);
  return false;
}

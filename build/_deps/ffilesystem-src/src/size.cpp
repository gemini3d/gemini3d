#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include "ffilesystem.h"

#include <string_view>
#include <system_error>

#include <cstdint>  // uintmax_t

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#else
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <dirent.h>  // opendir, readdir, closedir
#endif

#if __has_include(<fcntl.h>)
#include <fcntl.h>   // AT_* constants for statx
#endif

#endif


std::uintmax_t fs_file_size(std::string_view path)
{
  // fs_file_size() like std::filesystem::file_size() is only for files, not directories, which are considered to have no size.
  // Returns (uintmax_t)(-1) on error, and sets errno or std::error_code.
  // different platforms treat non-file's size differently, so we need the fs_is_file() for consistency.

  std::error_code ec;

  if(!fs_is_file(path))
    return fs_unknown_size;

#if defined(HAVE_CXX_FILESYSTEM)

  auto s = Filesystem::file_size(path, ec);
  if (ec)
    fs_print_error(path, ec);

  return s;

#else

  int r = 0;
  const std::string cpath(path);
#if defined(HAVE_STATX)
  struct statx sx;
  r = ::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT, STATX_SIZE, &sx);
  if (r == 0)
    return sx.stx_size;
#endif

  if (r == 0 || errno == ENOSYS){
    if (struct stat s; !::stat(cpath.c_str(), &s))
      return s.st_size;
  }

#endif

  fs_print_error(path, ec);
  return fs_unknown_size;
}


bool fs_is_empty(std::string_view path)
{
  // directory or file empty
  // returns false if path doesn't exist

  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)
  if (bool e = Filesystem::is_empty(path, ec); !ec)
    return e;
#else

  const std::string cpath(path);

  if (!fs_is_dir(path))
    return fs_is_file(path) &&fs_file_size(path) == 0;

  // directory empty
#if defined(_WIN32)
  // https://docs.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-findfirstfilea
  WIN32_FIND_DATAA ffd;
  HANDLE hFind = FindFirstFileA((cpath + "/*").c_str(), &ffd);
  if (hFind == INVALID_HANDLE_VALUE) {
    fs_print_error(path);
    return false;
  }

  // RAII for FindClose
  struct FindHandleCloser {
    HANDLE h;
    ~FindHandleCloser() { if (h != INVALID_HANDLE_VALUE) FindClose(h); }
  } _closer{hFind};

  do
  {
      if(fs_trace) std::cout << "TRACE: is_empty: do " << ffd.cFileName << "\n";

      if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
        // std::set is much slower than a simple if
        if (std::string_view n(ffd.cFileName); n == "." || n == "..")
          continue;

      // directory that is not . or ..
        return false;
      }
      // any non-directory
      return false;
  } while (FindNextFileA(hFind, &ffd));

  // empty directory
  return true;
#else
// https://www.man7.org/linux/man-pages/man3/opendir.3.html
// https://www.man7.org/linux/man-pages/man3/readdir.3.html
// https://www.man7.org/linux/man-pages/man3/closedir.3.html
// https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/readdir.3.html

  if (DIR *d = ::opendir(cpath.c_str()); d)
  {
    // RAII for closedir
    struct DirCloser { DIR* d; ~DirCloser(){ if(d) ::closedir(d); } } _dc{d};
    struct dirent *entry;
  while ((entry = ::readdir(d)))
  {
#ifdef _DIRENT_HAVE_D_TYPE
    if (entry->d_type == DT_DIR)
#else
    if (fs_is_dir(cpath + "/" + entry->d_name))
#endif
    {
      if (std::string_view n(entry->d_name); n == "." || n == "..")
        continue;
      // directory that is not . or ..
      return false;
    }
    // any non-directory
    return false;
  }
    // empty directory
    return true;
  }

#endif

#endif

  fs_print_error(path, ec);
  return false;

}

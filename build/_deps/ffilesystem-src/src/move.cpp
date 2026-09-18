#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#else

#include <cstdio> // for std::remove, std::rename

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h> // for DeleteFile
#endif
#endif

#include <string_view>
#include <system_error>         // for error_code

#include "ffilesystem.h"

#include <iostream>

bool
fs_remove(std::string_view path)
{
  // remove a file or empty directory
  // returns false if path does not exist or is a non-empty directory
  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)
  // https://en.cppreference.com/w/cpp/filesystem/remove
  if(Filesystem::remove(path, ec) && !ec) FFS_LIKELY
    return true;
#else
  const std::string cpath(path);
  // https://en.cppreference.com/w/cpp/io/c/remove
  if(std::remove(cpath.c_str()) == 0)
    return true;

#if defined(_WIN32)
  // Windows std::remove is confused by symlink directories, not deleting them, so give an _unlink
  if(fs_exists(path)){
  // https://docs.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-deletefilea
  // Need RemoveDirectoryW when deleting a directory symlink with std::remove()
  // https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-createsymboliclinka#remarks
    if(RemoveDirectoryW(fs_win32_to_wide(path).c_str()) != 0)
      return true;
  }
#endif

#endif

  fs_print_error(path, ec);
  return false;
}


bool
fs_rename(std::string_view from, std::string_view to)
{
  // rename a file or directory
  // existing "to" is overwritten

  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)
  // https://en.cppreference.com/w/cpp/filesystem/rename
  if(Filesystem::rename(from, to, ec); !ec)
#else
  // https://en.cppreference.com/w/cpp/io/c/rename
  const std::string cfrom(from), cto(to);
  if(std::rename(cfrom.c_str(), cto.c_str()) == 0)
#endif
    return true;

  fs_print_error(from, to, ec);
  return false;

}

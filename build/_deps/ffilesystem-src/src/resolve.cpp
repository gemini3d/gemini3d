#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE // for realpath
#endif
#endif

#include "ffilesystem.h"

#include <string>
#include <string_view>
#include <system_error>

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#endif

#if !defined(_WIN32)
#include <vector>
#endif


std::string
fs_canonical(
  std::string_view path,
  const bool strict)
{
  // canonicalize path, i.e. resolve all symbolic links, remove ".", ".." and extra slashes
  // if strict is true, then path must exist

  if (path.empty())
    return {};
    // need this for macOS otherwise it returns the current working directory instead of empty string

  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)

  if (fs_is_mingw() && fs_is_symlink(path))
    return fs_win32_final_path(path);

  if(auto c = strict ? Filesystem::canonical(path, ec) : Filesystem::weakly_canonical(path, ec); !ec)
    return c.generic_string();

#else

  if (std::string c = fs_realpath(path); !c.empty())
    return c;

  if (!strict)
    return fs_normal(path);
  // std::filesystem::weakly_canonical() and boost::filesystem::canonical()
  // resolve the path to the current working directory for non-existing paths.
  // https://github.com/boostorg/filesystem/blob/f4bb6d0f3ebe9f8b90243d8a98989191925d49d2/src/operations.cpp#L5023
  // That has significant complexity we will not implement here, as this is just a fallback method.
  // We just return the normal path.

#endif

  fs_print_error(path, ec);
  return {};
}


std::string fs_resolve(std::string_view path, const bool strict)
{
  // works like canonical(absolute(path)).
  // Inspired by Python pathlib.Path.resolve()
  // https://docs.python.org/3/library/pathlib.html#pathlib.Path.resolve

  return fs_canonical(fs_absolute(path), strict);
}


std::string fs_realpath(std::string_view path)
{
  // resolve real existing path
  // not defined for non-existing path--may return empty string.
  // https://man7.org/linux/man-pages/man3/realpath.3.html
  // https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/realpath.3.html

#if defined(_WIN32)
  return fs_win32_final_path(path);
#else
  std::vector<char> buf(fs_get_max_path());

  if (std::string cpath(path); ::realpath(cpath.c_str(), buf.data()))
    return buf.data();

  return {};
#endif

}

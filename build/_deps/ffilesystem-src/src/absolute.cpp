#include <string>
#include <string_view>

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
#include <system_error>
namespace Filesystem = std::filesystem;
#endif

#include "ffilesystem.h"


std::string fs_absolute(std::string_view path)
{
  // wraps std::filesystem::absolute(path).
  // path need not exist
  // Inspired by Python pathlib.Path.absolute()
  // https://docs.python.org/3/library/pathlib.html#pathlib.Path.absolute

  // Linux, MinGW can't handle empty paths
  if(path.empty())
    return fs_get_cwd();

  if (fs_is_absolute(path))
    return std::string(path);

#ifdef HAVE_CXX_FILESYSTEM
  std::error_code ec;

  if(auto a = Filesystem::absolute(path, ec); !ec)
    return a.string();

  fs_print_error(path, ec);
  return {};
#else
  std::string a = fs_get_cwd();
  if(a.empty())
    return {};

  a.push_back(fs_filesep());
  a += path;
  return a;
  // NOT normalized to be consistent with <filesystem>
#endif
}



std::string fs_absolute(std::string_view path, std::string_view base)
{
  // rebase path on base.

  if(fs_is_absolute(path))
    return std::string(path);

  std::string b = fs_absolute(base);
  if(b.empty())  FFS_UNLIKELY
    return {};

  if(!path.empty()){
    if(b.back() != '/' && b.back() != fs_filesep())
      b.push_back(fs_filesep());

    // don't need join(). Keeps it like Python pathlib.Path.absolute()
    b += path;
  }

  return b;
}

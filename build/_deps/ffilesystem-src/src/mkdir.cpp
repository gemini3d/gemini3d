#include "ffilesystem.h"

#include <string>  // IWYU pragma: keep
#include <string_view>

#include <iostream>
#include <system_error>         // for error_code
#include <vector>  // IWYU pragma: keep

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#else
#include <cerrno>
#include <sys/types.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#include <sys/stat.h>  // mkdir
#endif

#endif


bool fs_mkdir(std::string_view path)
{
 // make directory with parents

  std::error_code ec;

#ifdef HAVE_CXX_FILESYSTEM
  // https://en.cppreference.com/w/cpp/filesystem/create_directory
  Filesystem::create_directories(path, ec);
  // we don't use the return value as it indicates if a directory was created or not.
  if (!ec) FFS_LIKELY
    return true;
#else

if(path.empty())
  ec = std::make_error_code(std::errc::invalid_argument);
else {

  std::string buf;

  std::string const p = fs_resolve(path, false);
  // ERROR_PATH_NOT_FOUND if relative directory

  const std::vector<std::string> parts = fs_split(p);

  // if first part is root
  if(fs_slash_first(p))
    buf.push_back('/');

  // iterate over parts
  bool ok = false;

  for(const auto& p : parts){
    buf += p;
    buf.push_back('/');

    if(fs_trace) std::cout << "TRACE:mkdir " << buf << "\n";
    // create directory
#if defined(_WIN32)
    // https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-createdirectoryw
    ok = CreateDirectoryW(fs_win32_to_wide(buf).c_str(), nullptr) != 0;
    if(!ok){
      const auto err = GetLastError();
      // ERROR_PATH_NOT_FOUND if relative directory, thus we absolute() before loop.
      // ERROR_ACCESS_DENIED is OK if it's not the last part
      ok = err == ERROR_ALREADY_EXISTS ||
            (err == ERROR_ACCESS_DENIED && p != parts.back());
    }
#else
  // https://www.man7.org/linux/man-pages/man2/mkdir.2.html
    ok = ::mkdir(buf.c_str(), S_IRWXU) == 0 ||
           errno == EEXIST ||
          (errno == EACCES && p != parts.back());
#endif

    if(!ok)  FFS_UNLIKELY
      break;
  } // for

  if(ok)  FFS_LIKELY
    return true;

} // if path.empty()

#endif

  fs_print_error(path, ec);
  return false;
}

#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include "ffilesystem.h"

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#else
#include <vector>
#endif

#include <string>
#include <string_view>



#include <system_error> // for std::error_code

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h> // IWYU pragma: keep
// GetTempPathA
#else
#include <unistd.h> // IWYU pragma: keep
// getcwd, chdir
#endif


bool fs_set_cwd(std::string_view path)
{

  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)
  if(Filesystem::current_path(path, ec); !ec)
    return true;
#elif defined(_WIN32)
  // windows.h https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-setcurrentdirectory
  if(SetCurrentDirectoryW(fs_win32_to_wide(path).c_str()))
    return true;
#else
  // unistd.h https://www.man7.org/linux/man-pages/man2/chdir.2.html
  if(std::string cpath(path); ::chdir(cpath.c_str()) == 0)
    return true;
#endif

  fs_print_error(path, ec);
  return false;
}


std::string fs_get_cwd()
{
  // does not have trailing slash
  std::error_code ec;

#if defined(HAVE_CXX_FILESYSTEM)
  if(auto s = Filesystem::current_path(ec); !ec)
    return s.string();
#elif defined(_WIN32)
// windows.h https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-getcurrentdirectory
  if(DWORD L = GetCurrentDirectoryW(0, nullptr); L > 0) {
    std::wstring w;
    w.resize(L);
    if(GetCurrentDirectoryW(L, w.data()) == L-1) {
      w.resize(L-1);
      return fs_win32_to_narrow(w);
    }
  }
#else
// unistd.h https://www.man7.org/linux/man-pages/man3/getcwd.3.html
// https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/getcwd.3.html

std::vector<char> buf(fs_get_max_path());

if (::getcwd(buf.data(), buf.size()))
  return buf.data();

#endif

  fs_print_error("", ec);
  return {};
}

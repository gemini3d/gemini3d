#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE // for readlink
#endif
#endif

#include "ffilesystem.h"

#include <string>
#include <system_error>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h> // GetModuleFileName
#include <cstddef> // for size_t
#elif defined(FFS_DARWIN)
#include <cstdint> // for uint32_t
#include <mach-o/dyld.h> // _NSGetExecutablePath
#elif defined(__linux__) || defined(__CYGWIN__)
#include <unistd.h> // for readlink
#include <sys/types.h> // for ssize_t
#elif defined(FFS_BSD)
#include <sys/sysctl.h> // for sysctl
#include <cstddef> // for size_t
#endif


std::string fs_exe_path()
{
  // https://stackoverflow.com/a/4031835
  // https://stackoverflow.com/a/1024937

  std::error_code ec;

#if defined(_WIN32)
  // https://learn.microsoft.com/en-us/windows/win32/api/libloaderapi/nf-libloaderapi-getmodulefilenamea
  std::wstring w;
  w.resize(fs_get_max_path());

  if (DWORD L = GetModuleFileNameW(nullptr, w.data(), static_cast<DWORD>(w.size())); L > 0 && GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
    w.resize(L);
    return fs_win32_to_narrow(w);
  }
#elif defined(__linux__) || defined(__CYGWIN__)
  // https://man7.org/linux/man-pages/man2/readlink.2.html
  const std::string exe = "/proc/self/exe";
  std::string p;
  p.resize(fs_symlink_length(exe));

  if(ssize_t L = ::readlink(exe.c_str(), p.data(), p.size()); L > 0) {
    p.resize(L);
    return p;
  }
#elif defined(FFS_DARWIN)
  // https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/dyld.3.html
  // get buffer size first. Need L=0 to avoid intermittent segfault.
  std::uint32_t L = 0;

  if(_NSGetExecutablePath(nullptr, &L) == -1) {
    std::string path;
    path.resize(L);
    if(_NSGetExecutablePath(path.data(), &L) == 0){
      path.resize(L-1);
      return path;
    }
  }
#elif defined(FFS_BSD)
  // https://man.freebsd.org/cgi/man.cgi?sysctl(3)
  auto L = fs_get_max_path();
  std::string path;
  path.resize(L);

  const int mib[4] = {CTL_KERN, KERN_PROC, KERN_PROC_PATHNAME, -1};

  if(sysctl(mib, 4, path.data(), &L, nullptr, 0) == 0) {
    path.resize(L-1);
    return path;
  }
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error("", ec);
  return {};
}

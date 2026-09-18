#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#if defined(_WIN32)
#include <winsock2.h>
#elif __has_include(<sys/utsname.h>)
#define HAVE_UTSNAME
#include <sys/utsname.h>
#endif

#include <system_error>

#include <string>

#include "ffilesystem.h"


std::string fs_hostname()
{

  std::error_code ec;

#if defined(_WIN32)
  // https://learn.microsoft.com/en-us/windows/win32/api/winsock/nf-winsock-gethostname

  if(WSADATA wsaData; !WSAStartup(MAKEWORD(2, 0), &wsaData)){
    std::string name;
    name.resize(256);

    int r = gethostname(name.data(), static_cast<int>(name.size()));

    WSACleanup();

    if (r == 0) {
      fs_trim(name);
      return name;
    }
  }

#elif defined(HAVE_UTSNAME)
  // gethostname() doesn't dynamically allocate buffers, so use uname,
  // which gives the same result.

  if (struct utsname s; ::uname(&s) == 0)
    return s.nodename;
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error("", ec);
  return {};
}

#include <string>
#include <string_view>
#include <cstring> // for std::strcmp

#include <system_error>

#include "ffilesystem.h"

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#if __has_include(<format>)
#include <format> // IWYU pragma: keep
#endif

#endif

#if __has_include(<sys/utsname.h>)
#define HAVE_UTSNAME
#include <sys/utsname.h>
#endif


int fs_is_wsl()
{
  // return 0 if not WSL, 1 if WSL1, 2 if WSL2, -1 on error

#ifdef HAVE_UTSNAME
  struct utsname b;
  if (::uname(&b) != 0)
    return -1;

  if(std::strcmp(b.sysname, "Linux") != 0)
    return 0;

  std::string_view r(b.release);

#ifdef __cpp_lib_starts_ends_with // C++20
  if (r.ends_with("microsoft-standard-WSL2"))
    return 2;
  if (r.ends_with("-Microsoft"))
    return 1;
#endif
#endif
  return 0;
}


std::string fs_cpu_arch()
{

  std::error_code ec;

#ifdef HAVE_UTSNAME
  if (struct utsname b; ::uname(&b) == 0)
    return b.machine;
#elif defined(_WIN32)
// https://learn.microsoft.com/en-us/windows/win32/api/sysinfoapi/ns-sysinfoapi-system_info
    SYSTEM_INFO si;
    GetNativeSystemInfo(&si);
    switch (si.wProcessorArchitecture) {
    case PROCESSOR_ARCHITECTURE_AMD64: return "x86_64";
    case PROCESSOR_ARCHITECTURE_ARM64: return "arm64";
    case PROCESSOR_ARCHITECTURE_ARM: return "arm";
    case PROCESSOR_ARCHITECTURE_INTEL: return "x86";
    default: return "unknown";
    }
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error("", ec);
  return {};
}


std::string fs_os_version()
{
// get operating system version

  std::error_code ec;

#ifdef HAVE_UTSNAME
  if (struct utsname buf; ::uname(&buf) == 0)
    return buf.version;
#elif defined(_WIN32)

  OSVERSIONINFOA osvi;
  ZeroMemory(&osvi, sizeof(OSVERSIONINFOA));
  osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOA);

  if(GetVersionExA(&osvi))
#if defined(__cpp_lib_format)  // C++20
    return std::format("{}.{}.{}", osvi.dwMajorVersion, osvi.dwMinorVersion, osvi.dwBuildNumber);
#else
    return std::to_string(osvi.dwMajorVersion) + '.' + std::to_string(osvi.dwMinorVersion) + '.' + std::to_string(osvi.dwBuildNumber);
#endif

#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error("", ec);
  return {};
}

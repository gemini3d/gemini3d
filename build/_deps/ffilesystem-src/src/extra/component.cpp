#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include <string>
#include <string_view>

#include <system_error>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <limits.h>   // NAME_MAX
#include <unistd.h>  // pathconf
#include <cerrno>
#endif

#include "ffilesystem.h"

[[maybe_unused]] constexpr std::size_t DEFAULT_MAX_PATH = 256;

// This function returns the maximum length of a single component in the given path.
std::string_view::size_type fs_max_component(std::string_view path)
{
  // maximum length of each component of a path. That is, while the maximum
  // total length of a path may be thousands of character, each segment of the
  // path has a maximum length too--typically 255 characters.

  std::error_code ec = std::make_error_code(std::errc::invalid_argument);

#if defined(_WIN32)
  DWORD L = 0;
  if(GetVolumeInformationW(fs_win32_to_wide(fs_root(path)).c_str(), nullptr, 0, nullptr, &L, nullptr, nullptr, 0) != 0)
    return L;
#elif defined(_PC_NAME_MAX)
  errno = 0;
  const std::string cpath(path);
  auto const r = pathconf(cpath.c_str(), _PC_NAME_MAX);
  if(r != -1)
    return r;

  if(errno == 0)
    return DEFAULT_MAX_PATH;
#elif defined(NAME_MAX)
  return NAME_MAX;
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(path, ec);
  return {};
}

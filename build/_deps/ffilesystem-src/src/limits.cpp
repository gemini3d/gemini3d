#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include "ffilesystem.h"

#if defined(FFS_DARWIN)
#include <sys/syslimits.h>
#elif defined(_WIN32)
#include <cstdlib> // for _MAX_PATH
#endif

#if __has_include(<limits.h>)
#include <limits.h>
#endif



constexpr std::size_t DEFAULT_MAX_PATH = 256;

std::size_t
fs_get_max_path()
{
// Returns the maximum path length supported by the file system.
  auto m = DEFAULT_MAX_PATH;
#if defined(PATH_MAX)
  // POSIX
  m = PATH_MAX;
#elif defined(_MAX_PATH)
  // https://learn.microsoft.com/en-us/cpp/c-runtime-library/path-field-limits
  m = _MAX_PATH;
#elif defined(_POSIX_PATH_MAX)
  m = _POSIX_PATH_MAX;
#endif
  // arbitrary absolute maximum
  // Ref: https://github.com/gulrak/filesystem/blob/b1982f06c84f08a99fb90bac43c2d03712efe921/include/ghc/filesystem.hpp#L244
  return (m < 4096) ? m : 4096;
}

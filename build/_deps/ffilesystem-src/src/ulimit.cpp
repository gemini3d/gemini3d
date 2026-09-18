
#if defined(_WIN32)
#include <cstdio> // for _getmaxstdio
#else
#include <unistd.h>
#endif

#include "ffilesystem.h"


long fs_get_max_open_files()
// maximum number of open files
{
  long omax;
#if defined(_WIN32)
// https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/getmaxstdio?view=msvc-170
  omax = static_cast<long>(_getmaxstdio());
#else
  omax = ::sysconf(_SC_OPEN_MAX);
  if (omax < 0) {
    fs_print_error("", "sysconf");
  }
#endif
  return omax;
}

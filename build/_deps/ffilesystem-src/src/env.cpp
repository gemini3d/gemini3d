#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#include <cstdlib> // putenv, setenv, getenv

#include <string>
#include <string_view>
#include <optional>

#include "ffilesystem.h"


std::optional<std::string> fs_getenv(std::string_view name)
{
  // convenience function to get environment variable without needing to check for nullptr
  // don't emit error because sometimes we just check if envvar is defined
  const std::string n(name);
  auto buf = std::getenv(n.c_str());

  if (buf)
    return std::string(buf);

  return {};
}


bool fs_setenv(std::string_view name, std::string_view value)
{
  // if value is empty, remove the variable from the environment

const std::string n(name);
const std::string v(value);

#if defined(_WIN32)
  // SetEnvironmentVariable returned OK but set blank values
  // https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/putenv-wputenv
  if(_putenv_s(n.c_str(), v.c_str()) == 0)
    return true;
#else
  // https://www.man7.org/linux/man-pages/man3/setenv.3.html
  if(v.empty()) {
    if(::unsetenv(n.c_str()) == 0)
      return true;
  } else if(::setenv(n.c_str(), v.c_str(), 1) == 0)
    return true;
#endif

  fs_print_error(name);
  return false;
}

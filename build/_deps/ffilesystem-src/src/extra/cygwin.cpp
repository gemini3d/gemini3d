#include <string>
#include <string_view>

#include "ffilesystem.h"

#ifdef __CYGWIN__
#include <sys/cygwin.h>
#else
constexpr int CCP_WIN_A_TO_POSIX = 0;
constexpr int CCP_POSIX_TO_WIN_A = 0;
#endif


static std::string fs_convert_path(std::string_view path, [[maybe_unused]] int const what)
{
#ifdef __CYGWIN__
  const std::string cpath(path);
  const auto L = cygwin_conv_path(what, cpath.c_str(), nullptr, 0);
  if(L > 0){
    std::string r;
    r.resize(L);

    if (!cygwin_conv_path(what, cpath.c_str(), r.data(), L))
      return r;
  }
#endif

  fs_print_error(path);
  return {};
}

std::string fs_to_cygpath(std::string_view path) {
  return fs_convert_path(path, CCP_WIN_A_TO_POSIX);
}

std::string fs_to_winpath(std::string_view path) {
  return fs_convert_path(path, CCP_POSIX_TO_WIN_A);
}

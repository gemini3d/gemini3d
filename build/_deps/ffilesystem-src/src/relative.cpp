#include "ffilesystem.h"

#include <string>
#include <string_view>

#include <system_error>

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#endif


std::string fs_relative_to(std::string_view base, std::string_view other)
{
  // find relative path from base to other
  // "base" and "other" are treated as weakly canonical paths

  std::error_code ec;

  if (base.empty() && other.empty())
    return "";

#ifdef HAVE_CXX_FILESYSTEM
  if(std::string r = Filesystem::relative(other, base, ec).generic_string(); !ec)
    return r;
  // this implements
  // std::filesystem::weakly_canonical(other, ec).lexically_relative(std::filesystem::weakly_canonical(base, ec));
#else
   ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(base, other, ec);
  return {};
}


std::string fs_proximate_to(std::string_view base, std::string_view other)
{
  // find proximate path from base to other
  // "base" and "other" are treated as weakly canonical paths

  std::error_code ec;

  if (base.empty() && other.empty())
    return "";

#ifdef HAVE_CXX_FILESYSTEM
  if(std::string r = Filesystem::proximate(other, base, ec).generic_string(); !ec)
    return r;
#else
   ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(base, other, ec);
  return {};
}

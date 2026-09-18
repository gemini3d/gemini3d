#include <string>

#include "ffilesystem.h"

#if __has_include(<format>)
#include <format>
#endif
#if __has_include(<ranges>)
#include <ranges>
#endif


long fs_cpp_lang() {
  // C++ version compiler claims to support with given options
  return __cplusplus;
}


long fs_cpp_format() {
#ifdef __cpp_lib_format
  return __cpp_lib_format;
#else
  return 0;
#endif
}


long fs_cpp_ranges() {
#ifdef __cpp_lib_ranges
  return __cpp_lib_ranges;
#else
  return 0;
#endif
}


std::string fs_backend() {
#ifdef HAVE_CXX_FILESYSTEM
  return "<filesystem>";
#else
  return "C";
#endif
}


bool fs_is_optimized() {
// This is a heuristic, trusting the build system or user to set NDEBUG if optimized.
// The NDEBUG macro is typically defined when optimizations are enabled to disable debugging code.
// It is a standard way to indicate that the code should be optimized and not include debug information.
#if defined(NDEBUG)
  return true;
#else
  return false;
#endif
}

#include "ffilesystem.h"


long fs_c_lang(){
  // C version compiler claims to support with given options
#ifdef __STDC_VERSION__
  return __STDC_VERSION__;
#else  // MSVC if /std: not set
  return 0L;
#endif
}


bool fs_is_android(){
#ifdef __ANDROID__
  return true;
#else
  return false;
#endif
}

bool fs_is_macos(){
// we don't use the much more specific TargetConditionals.h because on macOS SDK updates sometimes including this header
// breaks compile until Homebrew / GCC updates.
#if defined(FFS_DARWIN)
  return true;
#else
  return false;
#endif
}

bool fs_is_linux() {
#ifdef __linux__
  return true;
#else
  return false;
#endif
}

bool fs_is_unix() {
#ifdef __unix__
  return true;
#else
  return fs_is_macos();
#endif
}

bool fs_is_bsd() {
#if defined(FFS_BSD)
  return true;
#else
  return false;
#endif
}

bool fs_is_windows() {
#ifdef _WIN32
  return true;
#else
  return false;
#endif
}

bool fs_is_cygwin(){
#ifdef __CYGWIN__
  return true;
#else
  return false;
#endif
}

bool fs_is_mingw(){
#ifdef __MINGW32__
  return true;
#else
  return false;
#endif
}

bool fs_is_appleclang(){
#if defined(__clang__) && defined(__apple_build_version__)
  return true;
#else
  return false;
#endif
}

bool fs_is_clangcl(){
#if defined(__clang__) && defined(_MSC_VER)
  return true;
#else
  return false;
#endif
}

bool fs_is_msvc(){
#ifdef _MSC_VER
  return true;
#else
  return false;
#endif
}

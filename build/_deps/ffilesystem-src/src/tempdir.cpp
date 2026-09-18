#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#endif

#include "ffilesystem.h"

#include <string>
#include <string_view>

#include <system_error> // for std::error_code

#include <iostream> // IWYU pragma: keep

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h> // IWYU pragma: keep
// GetTempPathA
#elif defined(FFS_DARWIN)
#include <unistd.h> // for confstr
#endif



std::string fs_get_tempdir()
{
  std::error_code ec;

#ifdef HAVE_CXX_FILESYSTEM
  if(auto p = Filesystem::temp_directory_path(ec); !ec && !p.empty())
    return p.string();
#endif

#if defined(_WIN32)
  // GetTempPath2 is not in MSYS2. libuv etc. use GetTempPathW
  if(DWORD L = GetTempPathW(0, nullptr); L > 0) {
    std::wstring w;
    w.resize(L + 1);
    if(GetTempPathW(L, w.data()) == L) {
      w.resize(L);
      return fs_win32_to_narrow(w);
    }
    if (L > 0) {
      w.resize(L);
      return fs_win32_to_narrow(w);
    }
  }
#elif defined(FFS_DARWIN)
// https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/confstr.3.html
  size_t len = ::confstr(_CS_DARWIN_USER_TEMP_DIR, nullptr, 0);
  if (len > 1) {
    std::string t;
    t.resize(len);
    if(::confstr(_CS_DARWIN_USER_TEMP_DIR, t.data(), len) == len) {
      t.resize(len - 1); // remove trailing null
      if(fs_trace) std::cout << "TRACE: used confstr(_CS_DARWIN_USER_TEMP_DIR) = " << t << "\n";
      return t;
    }
  }
#endif

#if !defined(_WIN32)
  if(auto t = fs_getenv("TMPDIR"); t.has_value() && !t.value().empty())
    return t.value();

  if (fs_is_dir("/tmp"))
    return "/tmp";
#endif

  fs_print_error("", ec);
  return {};

}

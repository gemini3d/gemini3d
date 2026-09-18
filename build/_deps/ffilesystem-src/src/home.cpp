#include <string>
#include <string_view>

#include <system_error>

#if __has_include(<format>)
#include <format>  // IWYU pragma: keep
#endif


// get_profile_dir
#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#ifndef SECURITY_WIN32
#define SECURITY_WIN32
#endif
#include <UserEnv.h> // GetUserProfileDirectory
#include <Security.h> // GetUserNameEx
#include <Windows.h>
#else
#include <sys/types.h>  // IWYU pragma: keep
#include <pwd.h>      // for getpwuid, passwd
#include <unistd.h> // for mac too
#endif

#include "ffilesystem.h"


struct passwd* fs_getpwuid()
{
#if !defined(_WIN32)
  const uid_t eff_uid = ::geteuid();

  if(auto pw = ::getpwuid(eff_uid))
    return pw;

  fs_print_error(
#if defined(__cpp_lib_format)  // C++20
    std::format("uid: {}", eff_uid)
#else
    std::to_string(eff_uid)
#endif
    );
#endif

  return {};
}


std::string fs_get_homedir()
{
  if (auto h = fs_getenv(fs_is_windows() ? "USERPROFILE" : "HOME"); h.has_value() && !h.value().empty())
    return h.value();

  return fs_get_profile_dir();
}


std::string fs_get_profile_dir()
{

#if defined(_WIN32)
  // https://learn.microsoft.com/en-us/windows/win32/api/userenv/nf-userenv-getuserprofiledirectorya
  // works on MSYS2, MSVC, oneAPI
  HANDLE h = nullptr;

  if(OpenProcessToken( GetCurrentProcess(), TOKEN_QUERY, &h)) {
    DWORD L = 0;
    GetUserProfileDirectoryW(h, nullptr, &L);
    BOOL ok = false;
    std::wstring w;

    if(L > 0) {
      w.resize(L);
      ok = GetUserProfileDirectoryW(h, w.data(), &L);
    }

    CloseHandle(h);

    if(ok && L > 0)
      return fs_win32_to_narrow(w);
  }
#else
  if (auto pw = fs_getpwuid())
    return pw->pw_dir;
#endif

  fs_print_error("");
  return {};
}


std::string fs_expanduser(std::string_view path)
{
  if(path.empty())
    return {};

  if(path.front() != '~')
    return std::string(path);

  // second character is not a file separator
  // std::set is much slower than a simple if
  if(path.length() > 1 && !(path[1] == '/' || path[1] == fs_filesep()))
    return std::string(path);

  const std::string home = fs_get_homedir();
  if(home.empty())
    return {};

  if (path.length() < 3)
    return home;

// handle initial duplicated file separators. NOT .lexical_normal to handle "~/.."
// std::set is much slower than a simple if
  std::string::size_type i = 2;
  while(i < path.length() && (path[i] == '/' || path[i] == fs_filesep()))
    i++;

  std::string e = home;

  // e.reserve(home.size() + 1 + path.size() - i);
  // the .reserve makes no measureable improvement even at nanosecond scale.
  if (e.back() != '/' && e.back() != fs_filesep())
    e.push_back('/');

  e += path.substr(i);

  return e;
}


std::string fs_get_username()
{
  // Get username of the current user

#if defined(_WIN32)

// https://learn.microsoft.com/en-us/windows/win32/api/secext/nf-secext-getusernameexa
// https://learn.microsoft.com/en-us/windows/win32/api/secext/ne-secext-extended_name_format
  ULONG L = 0;
  if (GetUserNameExW(NameSamCompatible, nullptr, &L) == 0 && L > 0) {
    std::wstring w;
    w.resize(L);
    if (GetUserNameExW(NameSamCompatible, w.data(), &L) != 0)
      return fs_win32_to_narrow(w);
  }

#else

  if (auto pw = fs_getpwuid())
    return pw->pw_name;

#endif

  fs_print_error("");
  return {};
}

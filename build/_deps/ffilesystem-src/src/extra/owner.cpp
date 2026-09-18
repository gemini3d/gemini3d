// IWYU pragma: no_include <bits/statx-generic.h>
// IWYU pragma: no_include <linux/stat.h>

#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <aclapi.h> // GetNamedSecurityInfo
#include <Windows.h>
#else
#include <sys/types.h>  // IWYU pragma: keep
#include <sys/stat.h>  // IWYU pragma: keep
#include <pwd.h>      // for getpwuid
#include <grp.h>     // for getgrgid
#include <optional>
#endif

#if __has_include(<fcntl.h>)
#include <fcntl.h>   // AT_* constants for statx
#endif

#include <string>
#include <string_view>

#include "ffilesystem.h"


#if defined(_WIN32)

static std::string fs_win32_get_owner(PSID pSid)
{
  DWORD L1 = 0;
  DWORD L2 = 0;
  SID_NAME_USE eUse = SidTypeUnknown;

  if(!LookupAccountSidA(nullptr, pSid, nullptr, &L1, nullptr, &L2, &eUse) && GetLastError() != ERROR_INSUFFICIENT_BUFFER)
    return {};

  std::string s;
  s.resize(L1);

  if (!LookupAccountSidA(nullptr, pSid, s.data(), &L1, nullptr, &L2, &eUse))
    return {};

  // it's L1, not L1 - 1
  s.resize(L1);
  return s;
}

static std::string fs_win32_owner(std::string_view path, bool group)
{
// https://learn.microsoft.com/en-us/windows/win32/secauthz/finding-the-owner-of-a-file-object-in-c--

  PSECURITY_DESCRIPTOR pSD = nullptr;
  PSID pSid = nullptr;
  DWORD r;
  const std::string cpath(path);
  // https://learn.microsoft.com/en-us/windows/win32/api/aclapi/nf-aclapi-getnamedsecurityinfoa
  if (group)
    r = GetNamedSecurityInfoA(cpath.c_str(), SE_FILE_OBJECT, GROUP_SECURITY_INFORMATION, nullptr, &pSid, nullptr, nullptr, &pSD);
  else
    r = GetNamedSecurityInfoA(cpath.c_str(), SE_FILE_OBJECT, OWNER_SECURITY_INFORMATION, &pSid, nullptr, nullptr, nullptr, &pSD);

  std::string s;
  if(r == ERROR_SUCCESS)
    s = fs_win32_get_owner(pSid);

  LocalFree(pSD);
  // https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-localfree
  return s;
}


#else

static std::optional<uid_t> fs_stat_uid(std::string_view path)
{
  int r = 0;
  const std::string cpath(path);

#if defined(HAVE_STATX)
  struct statx sx;
  r = ::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT, STATX_UID, &sx);
  if (r == 0)
    return sx.stx_uid;
#endif

  if(r == 0 || errno == ENOSYS){
    if(struct stat s; !::stat(cpath.c_str(), &s))
      return s.st_uid;
  }

  return {};
}

static std::optional<gid_t> fs_stat_gid(std::string_view path)
{
  int r = 0;
  const std::string cpath(path);

#if defined(HAVE_STATX)
  struct statx sx;
  r = ::statx(AT_FDCWD, cpath.c_str(), AT_NO_AUTOMOUNT, STATX_GID, &sx);
  if (r == 0)
    return sx.stx_gid;
#endif

if(r == 0 || errno == ENOSYS){
  if(struct stat s; !::stat(cpath.c_str(), &s))
    return s.st_gid;
}

  return {};
}

#endif


std::string
fs_get_owner_name(std::string_view path)
{
#if defined(_WIN32)
  if (std::string s = fs_win32_owner(path, false); !s.empty())
    return s;
#else
  if (auto uid = fs_stat_uid(path)) {
    if (auto pw = getpwuid(uid.value()))
      return pw->pw_name;
  }
#endif

  fs_print_error(path);
  return {};
}


std::string
fs_get_owner_group(std::string_view path)
{
#if defined(_WIN32)
  if (std::string s = fs_win32_owner(path, true); !s.empty())
    return s;
#else
  if (auto gid = fs_stat_gid(path)) {
    if (auto gr = getgrgid(gid.value()))
      return gr->gr_name;
  }
#endif

  fs_print_error(path);
  return {};
}

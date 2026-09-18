#if defined(_WIN32)
# define WIN32_LEAN_AND_MEAN
# include <windows.h>  // GetTokenInformation
# include <io.h> // _isatty
#else
#if defined(__linux__) || defined(__CYGWIN__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#endif
# include <unistd.h>  // geteuid, getpid, isatty
# include <sys/types.h>  // IWYU pragma: keep
// geteuid, pid_t
#endif

#include <cstdio> // fileno

#include <string>

#include "ffilesystem.h"


bool fs_is_admin(){
  // is running as admin / root / superuser
#if defined(_WIN32)
	HANDLE h = nullptr;
	TOKEN_ELEVATION elevation;
	DWORD dwSize;

  // https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-openprocesstoken
  // https://learn.microsoft.com/en-us/windows/win32/api/securitybaseapi/nf-securitybaseapi-gettokeninformation

  const bool ok = (OpenProcessToken(GetCurrentProcess(), TOKEN_QUERY, &h) &&
     GetTokenInformation(h, TokenElevation, &elevation, sizeof(elevation), &dwSize));

  if (h)
    CloseHandle(h);

  if(ok)
    return elevation.TokenIsElevated;

  fs_print_error("");
  return false;
#else
  return ::geteuid() == 0;
#endif
}


pid_t fs_getpid()
{
  // get process ID
#if defined(_WIN32)
  // https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-getcurrentprocessid
  return static_cast<pid_t>(GetCurrentProcessId());
#else
  // https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man2/getpid.2.html
  // https://www.man7.org/linux/man-pages/man2/getpid.2.html
  return ::getpid();
#endif
}


bool fs_stdin_tty()
{
// detect if stdin is connected to a terminal
  return
#if defined(_WIN32)
    ::_isatty(::_fileno(stdin));
#else
    ::isatty(::fileno(stdin));
#endif
}


std::string fs_get_terminal()
{
#if defined(_WIN32)
  // inspired by https://gitlab.kitware.com/utils/kwsys/-/commit/0d6eac1feb8615fe59e8f972d41d1eaa8bc9baf8
  // Windows Console Host: ConsoleWindowClass
  // Windows Terminal / ConPTY: PseudoConsoleWindow (undocumented)

  // https://learn.microsoft.com/en-us/windows/console/getconsolewindow
  // encourages Virtual Terminal Sequences

  // The `GetConsoleWindow()` API returns a window handle (HWND), but it should not be closed by the application.
  // This handle is owned and managed by the system, so calling CloseWindow() on it could
  // disrupt the console window's operation.

  std::string name;
  name.resize(fs_get_max_path());

  if (HWND h = GetConsoleWindow(); h) {
    // https://learn.microsoft.com/en-us/windows/win32/api/winuser/nf-winuser-getclassnamea
    int L = GetClassNameA(h, name.data(), static_cast<int>(name.size()));

    if (L > 0){
      name.resize(L);
      return name;
    }
  }

#else
  if (auto t = fs_getenv("TERM"); t)
    return t.value();
#endif

  fs_print_error("");
  return {};
}

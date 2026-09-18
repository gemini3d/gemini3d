#if defined(_WIN32) || defined(__CYGWIN__)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winioctl.h>
#endif

#include <cstddef>  // for std::byte
#include <iostream>

#include <string>
#include <string_view>

#include <system_error>

#include "ffilesystem.h"

#if defined(_WIN32) || defined(__CYGWIN__)
// create type PREPARSE_DATA_BUFFER
// from ntifs.h, which can only be used by drivers
// typedef is copied from https://gitlab.kitware.com/utils/kwsys/-/blob/master/SystemTools.cxx
// that has a BSD 3-clause license
typedef struct _REPARSE_DATA_BUFFER
{
  ULONG ReparseTag;
  USHORT ReparseDataLength;
  USHORT Reserved;
  union
  {
    struct
    {
      USHORT SubstituteNameOffset;
      USHORT SubstituteNameLength;
      USHORT PrintNameOffset;
      USHORT PrintNameLength;
      ULONG Flags;
      WCHAR PathBuffer[1];
    } SymbolicLinkReparseBuffer;
    struct
    {
      USHORT SubstituteNameOffset;
      USHORT SubstituteNameLength;
      USHORT PrintNameOffset;
      USHORT PrintNameLength;
      WCHAR PathBuffer[1];
    } MountPointReparseBuffer;
    struct
    {
      UCHAR DataBuffer[1];
    } GenericReparseBuffer;
    struct
    {
      ULONG Version;
      WCHAR StringList[1];
      // In version 3, there are 4 NUL-terminated strings:
      // * Package ID
      // * Entry Point
      // * Executable Path
      // * Application Type
    } AppExecLinkReparseBuffer;
  } DUMMYUNIONNAME;
} REPARSE_DATA_BUFFER, *PREPARSE_DATA_BUFFER;
#endif


#if defined(_WIN32) || defined(__CYGWIN__)
static bool fs_win32_get_reparse_buffer(std::string_view path, std::byte* buffer)
{

// this function is adapted from
// https://gitlab.kitware.com/utils/kwsys/-/blob/master/SystemTools.cxx
// that has a BSD 3-clause license

  std::error_code ec;


// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getfileattributesw

// https://learn.microsoft.com/en-us/windows/win32/fileio/file-attribute-constants
  if (DWORD attr = GetFileAttributesW(fs_win32_to_wide(path).c_str()); attr == INVALID_FILE_ATTRIBUTES)
    // ec = std::make_error_code(std::errc::no_such_file_or_directory);
    // don't emit error for non-existent files
    return false;
  else if (!(attr & FILE_ATTRIBUTE_REPARSE_POINT))
    return false;

  // Using 0 instead of GENERIC_READ as it allows reading of file attributes
  // even if we do not have permission to read the file itself

  // A reparse point may be an execution alias (Windows Store app), which
  // is similar to a symlink but it cannot be opened as a regular file.
  // We must look at the reparse point data explicitly.

  // FILE_ATTRIBUTE_REPARSE_POINT means:
  // * a file or directory that has an associated reparse point, or
  // * a file that is a symbolic link.

  HANDLE h = CreateFileW(fs_win32_to_wide(path).c_str(),
                         0, 0, nullptr, OPEN_EXISTING,
                         FILE_FLAG_OPEN_REPARSE_POINT | FILE_FLAG_BACKUP_SEMANTICS, nullptr);

  if (h == INVALID_HANDLE_VALUE)
    ec = std::make_error_code(std::errc::io_error);
  else {

    DWORD bytesReturned = 0;

    BOOL ok = DeviceIoControl(h, FSCTL_GET_REPARSE_POINT, nullptr, 0, buffer,
                          MAXIMUM_REPARSE_DATA_BUFFER_SIZE, &bytesReturned,
                          nullptr);

    CloseHandle(h);
    if(ok)
      return true;
  }


  fs_print_error(path, ec);
  return false;
}
#endif


bool fs_is_appexec_alias(std::string_view path)
{
// Windows App Execution Alias allow easy access to Windows Store apps
// from the Windows Terminal. However, they don't work with _access_s()
// or stat() like regular files. This function detects the aliases.
// Reference:
// https://learn.microsoft.com/en-us/windows/terminal/command-line-arguments?tabs=windows#add-windows-terminal-executable-to-your-path
// https://learn.microsoft.com/en-us/windows/win32/fileio/reparse-points
//
// this function is adapted from
// https://gitlab.kitware.com/utils/kwsys/-/blob/master/SystemTools.cxx
// that has a BSD 3-clause license

#if defined(_WIN32) || defined(__CYGWIN__)
  std::byte buffer[MAXIMUM_REPARSE_DATA_BUFFER_SIZE];

  if (!fs_win32_get_reparse_buffer(path, buffer))
    return false;

  PREPARSE_DATA_BUFFER data =
    reinterpret_cast<PREPARSE_DATA_BUFFER>(&buffer[0]);

  // https://learn.microsoft.com/en-us/windows/win32/fileio/reparse-point-tags
  return data->ReparseTag == IO_REPARSE_TAG_APPEXECLINK;

#else
  fs_print_error(path, std::make_error_code(std::errc::function_not_supported));
  return false;
#endif

}


bool fs_win32_long_paths_enabled() {
  // from https://github.com/microsoft/STL/pull/5783/
  // microsoft/STL has Apache 2.0 license
  // https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation?tabs=powershell
  // https://learn.microsoft.com/en-us/windows/win32/api/winreg/nf-winreg-reggetvaluew

#if defined(_WIN32)
  DWORD registry_value = 0;
  DWORD buffer_size    = sizeof(registry_value);
  const LSTATUS status = RegGetValueW(HKEY_LOCAL_MACHINE, LR"(SYSTEM\CurrentControlSet\Control\FileSystem)",
      L"LongPathsEnabled", RRF_RT_REG_DWORD, nullptr, &registry_value, &buffer_size);

  switch (status) {
  case ERROR_SUCCESS:
      if(buffer_size == sizeof(registry_value))
        return registry_value != 0;
      break;
  case ERROR_FILE_NOT_FOUND:
      return false; // The registry value doesn't exist, so long paths aren't enabled.
  case ERROR_MORE_DATA:
      std::cerr << "fs_win32_long_paths_enabled: RegGetValueW() returned ERROR_MORE_DATA; this should not be possible.\n";
      break;
  case ERROR_UNSUPPORTED_TYPE:
      std::cerr << "fs_win32_long_paths_enabled: RegGetValueW() returned ERROR_UNSUPPORTED_TYPE; the value exists but has a weird type.\n";
      break;
  default:
      std::cerr << "fs_win32_long_paths_enabled: RegGetValueW() returned " << status << ".\n";
      break;
  }

#else
  fs_print_error("", std::make_error_code(std::errc::function_not_supported));
#endif

  return false;
}


bool fs_win32_is_symlink(std::string_view path)
{
// distinguish between Windows symbolic links and reparse points as
// reparse points can be unlike symlinks.
//
// this function is adapted from
// https://gitlab.kitware.com/utils/kwsys/-/blob/master/SystemTools.cxx
// that has a BSD 3-clause license

#if defined(_WIN32) || defined(__CYGWIN__)
  std::byte buffer[MAXIMUM_REPARSE_DATA_BUFFER_SIZE];
  // Since FILE_ATTRIBUTE_REPARSE_POINT is set this file must be
  // a symbolic link if it is not a reparse point.
  if (!fs_win32_get_reparse_buffer(path, buffer))
    return GetLastError() == ERROR_NOT_A_REPARSE_POINT;

  ULONG reparseTag =
    reinterpret_cast<PREPARSE_DATA_BUFFER>(&buffer[0])->ReparseTag;

  return (reparseTag == IO_REPARSE_TAG_SYMLINK) ||
         (reparseTag == IO_REPARSE_TAG_MOUNT_POINT);
#else
  fs_print_error(path, std::make_error_code(std::errc::function_not_supported));
  return false;
#endif

}


bool
fs_win32_is_ext_path(std::string_view path)
{
  // Windows extended path prefix enables paths longer than MAX_PATH (260 characters)
  // and disables path parsing, so that special characters are allowed

  //   \\?\ - enables extended-length paths and disables path parsing
constexpr std::string_view ext = R"(\\?\)";
  //   \\.\ - used for device paths (e.g., physical drives, COM ports)
constexpr std::string_view dev = R"(\\.\)";

  #ifdef __cpp_lib_starts_ends_with  // C++20
    return path.starts_with(ext) || path.starts_with(dev);
  #else
    if (path.length() < 4)
      return false;

    std::string_view p = path.substr(0, 4);
    return p == ext || p == dev;
  #endif
}


std::string fs_win32_full_name(std::string_view path)
{
// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getfullpathnameW
// GetFinalPathNameByHandle or GetFullPathName returns the unresolved symlink

  std::error_code ec;

#if defined(_WIN32) || defined(__CYGWIN__)
  std::wstring const w = fs_win32_to_wide(path);

  auto const L = GetFullPathNameW(w.c_str(), 0, nullptr, nullptr);
  // this form includes the null terminator
  // weak detection of race condition (cwd change)
  if(L){
    std::wstring r;
    r.resize(L);
    if(GetFullPathNameW(w.c_str(), L, r.data(), nullptr) == L-1)  FFS_LIKELY
    {
      r.resize(L-1);
      return fs_win32_to_narrow(r);
    }
  }
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(path, ec);
  return {};
}


std::string fs_win32_final_path(std::string_view path)
{
  // resolves Windows symbolic links (reparse points and junctions)
  // it also resolves the case insensitivity of Windows paths to the disk case
  // PATH MUST EXIST
  //
  // References:
  // https://stackoverflow.com/a/50182947
  // https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-createfilew
  // https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getfinalpathnamebyhandlew

  std::error_code ec;

#if defined(_WIN32) || defined(__CYGWIN__)
  if(fs_trace) std::cout << "TRACE: win32_final_path(" << path << ")\n";
  // dwDesiredAccess=0 to allow getting parameters even without read permission
  // FILE_FLAG_BACKUP_SEMANTICS required to open a directory

  std::wstring w = fs_win32_to_wide(path);

  HANDLE h = CreateFileW(w.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                         OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS, nullptr);
  if(h == INVALID_HANDLE_VALUE)
    return {};

  if(DWORD L = GetFinalPathNameByHandleW(h, nullptr, 0, FILE_NAME_NORMALIZED | VOLUME_NAME_DOS); L > 0) {
    w.resize(L + 1);

    L = GetFinalPathNameByHandleW(h, w.data(), L, FILE_NAME_NORMALIZED | VOLUME_NAME_DOS);

    CloseHandle(h);

    if(L > 0) {
      w.resize(L);
      std::string r(fs_win32_to_narrow(w));

      if (fs_win32_is_ext_path(r) && !fs_win32_is_ext_path(path))
        r = r.substr(4);

      return r;
    }
  } else {
    CloseHandle(h);
  }
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(path, ec);
  return {};
}


std::string fs_longname(std::string_view in)
{
// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getlongpathnamew
// the path must exist

  std::error_code ec;

#if defined(_WIN32) || defined(__CYGWIN__)
// size includes null terminator on 1st call, but does not include null terminator on 2nd call
  std::wstring const w = fs_win32_to_wide(in);
  DWORD L = GetLongPathNameW(w.c_str(), nullptr, 0);

  if(L > 0){
    std::wstring out;
    out.resize(L);

    if(GetLongPathNameW(w.c_str(), out.data(), L) == L-1) {
      out.resize(L);
      return fs_win32_to_narrow(out);
    }
  }
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(in, ec);
  return {};
}


std::string fs_shortname(std::string_view in)
{
// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getshortpathnamew
// the path must exist

  std::error_code ec;

#if defined(_WIN32) || defined(__CYGWIN__)
// size includes null terminator on 1st call, but does not include null terminator on 2nd call
  std::wstring const w = fs_win32_to_wide(in);
  DWORD L = GetShortPathNameW(w.c_str(), nullptr, 0);

  if(L > 0){
    std::wstring out;
    out.resize(L);

    if(GetShortPathNameW(w.c_str(), out.data(), L) == L-1) {
      out.resize(L);
      return fs_win32_to_narrow(out);
    }
  }
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(in, ec);
  return {};
}


std::string fs_win32_to_narrow([[maybe_unused]] std::wstring_view w)
{
  std::error_code ec;

#if defined(_WIN32) || defined(__CYGWIN__)
  std::wstring ws(w);
  if (int L = WideCharToMultiByte(CP_UTF8, 0, ws.c_str(), -1, nullptr, 0, nullptr, nullptr); L > 0)  FFS_LIKELY
  {
    std::string n;
    n.resize(L);

    if(WideCharToMultiByte(CP_UTF8, 0, ws.c_str(), -1, n.data(), L, nullptr, nullptr) == L) {
      n.resize(L-1);  // discard null terminator
      return n;
    }
  }
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error("", ec);
  return {};
}


std::wstring fs_win32_to_wide(std::string_view n)
{
  std::error_code ec;

#if defined(_WIN32) || defined(__CYGWIN__)
  // https://docs.microsoft.com/en-us/windows/win32/api/stringapiset/nf-stringapiset-multibytetowidechar
  std::string ns(n);
  if (int L = MultiByteToWideChar(CP_UTF8, 0, ns.c_str(), -1, nullptr, 0); L > 0)  FFS_LIKELY
  {
    std::wstring w;
    w.resize(L);

    if(MultiByteToWideChar(CP_UTF8, 0, ns.c_str(), -1, w.data(), L) == L) {
      w.resize(L-1);  // discard null terminator
      return w;
    }
  }
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(n, ec);
  return {};
}

#include <string>
#include <system_error>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <iostream>
#include <windows.h> // for OpenProcess, GetModuleBaseNameA
#include <tlhelp32.h> // for CreateToolhelp32Snapshot, Process32First
#include <psapi.h>  // for EnumProcessModules
#else
#include <pwd.h>     // for passwd
#endif

#include "ffilesystem.h"


std::string
fs_get_shell()
{

  std::error_code ec;

#if defined(_WIN32)
  const HANDLE h = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);
  if (h == INVALID_HANDLE_VALUE) {
    fs_print_error("");
    return {};
  }

  // 0: current process
  PROCESSENTRY32 pe;
  ZeroMemory(&pe, sizeof(PROCESSENTRY32));
  pe.dwSize = sizeof(PROCESSENTRY32);
  HMODULE hMod;
  DWORD cbNeeded;

  std::string name;
  name.resize(fs_get_max_path());

  if( Process32First(h, &pe)) {
    const DWORD pid = GetCurrentProcessId();
    do {
      if (pe.th32ProcessID == pid) {
        HANDLE hProcess = OpenProcess( PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, pe.th32ParentProcessID);
        if (!hProcess) {
          if (fs_trace) std::cout << "TRACE: OpenProcess failed for PPID " << pe.th32ParentProcessID << "\n";
          continue;
        }

        if ( EnumProcessModules( hProcess, &hMod, sizeof(hMod), &cbNeeded) ) {
          if(DWORD L = GetModuleBaseNameA( hProcess, hMod, name.data(), static_cast<DWORD>(name.size()) ); L == 0)
            ec = std::make_error_code(std::errc::io_error);
          else
            name.resize(L);

          CloseHandle(hProcess);
          break;
        }

        CloseHandle(hProcess);

        if(fs_trace) std::cout << "TRACE: get_shell: " << name << " PID: " << pid << " PPID: " << pe.th32ParentProcessID << "\n";
      }
    } while( Process32Next(h, &pe));
  }

  CloseHandle(h);

  if (!ec)
    return name;
#else
  if (auto pw = fs_getpwuid())
    return pw->pw_shell;
#endif

  fs_print_error("", ec);
  return {};
}

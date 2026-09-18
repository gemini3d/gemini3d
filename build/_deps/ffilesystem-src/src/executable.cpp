#include <string_view>
#include <string>
#include <cstdint> // for uint8_t

#include <iostream>

#include "ffilesystem.h"

// include even if <filesystem> is available
#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#else
#include <unistd.h>  // X_OK, access()
#endif

#if !defined(_WIN32)
#include <fstream>
#include <array>
#endif

#if defined(FFS_DARWIN)
#include <mach-o/loader.h>
#include <mach-o/fat.h>

#if __has_include(<bit>)
#include <bit>
#endif

#include <cstring>   // for std::memcpy. always include this for macOS to simplify logic and avoid build failures

#endif



bool fs_is_executable_binary(std::string_view path)
{
  // is path a binary executable file as determined by the magic number
  // beginning the file.  This is a heuristic and may not be 100% accurate.

  bool ok = false;
#if defined(_WIN32)
  // MinGW, MSVC, oneAPI at least need || is_appexec_alias()
  DWORD t;
  ok = (GetBinaryTypeW(fs_win32_to_wide(path).c_str(), &t) != 0) || fs_is_appexec_alias(path);
#elif defined(__CYGWIN__)
  fs_print_error(path, "not supported on Cygwin");
#else
  // https://github.com/jart/cosmopolitan/blob/master/ape/specification.md
  std::array<std::uint8_t, 4> magic;
  const std::string cpath(path);

  if(std::ifstream f{cpath.c_str(), std::ios::binary}){
    if( !f.read(reinterpret_cast<char*>(magic.data()), magic.size()) ) {
      return false;
    }
  } else {
    fs_print_error(path, "could not open file");
    return false;
  }

#if defined(FFS_DARWIN)
  // https://ss64.com/mac/file.html
    // Mach-O magic numbers
    // https://stackoverflow.com/q/27669766
    // https://jonasdevlieghere.com/post/macho-4gb-limit/
    // MH_MAGIC: 0xfeedface (32-bit)
    // MH_MAGIC_64: 0xfeedfacf (64-bit)
    // FAT_MAGIC: 0xCAFEBABE (32-bit)
    // FAT_MAGIC_64: 0xCAFEBABF (64-bit)
    // it's necessary to consider the reversed magic number for each of these as well,
    // done in the "*CIGAM*" constants
#if defined(__cpp_lib_bit_cast) // C++20
  auto um32 = std::bit_cast<std::uint32_t>(magic);
#else
  std::uint32_t um32;
  std::memcpy(&um32, magic.data(), sizeof(um32));
#endif

  ok = um32 == MH_MAGIC ||  um32 == MH_MAGIC_64 ||
       um32 == MH_CIGAM ||  um32 == MH_CIGAM_64 ||
       um32 == FAT_MAGIC || um32 == FAT_MAGIC_64 ||
       um32 == FAT_CIGAM || um32 == FAT_CIGAM_64;

#else
  // Linux / BSD: executable binary check via ELF file magic number
  // not using to_array to keep compatibility with C++17
  constexpr std::array<std::uint8_t, 4> ELF_MAGIC{0x7f, 'E', 'L', 'F'};

  ok = magic == ELF_MAGIC;
  // does not consider PE or COFF formats.
#endif
#endif

  return ok;
}


bool fs_is_exe(std::string_view path)
{
  // is file or symlink to a file executable by the user.
  // directories are not considered "executable file"
  // use fs_get_permissions() for that.

  if(!fs_is_file(path))
    return false;

#if defined(_WIN32)
  std::string suffix = fs_suffix(path);
  if (suffix.empty())
    return fs_is_executable_binary(path);

  if(!fs_is_readable(path))
    return false;

  std::string pathext = ".com;.exe;.bat;.cmd";
  if(auto e = fs_getenv("PATHEXT"); e) {
    pathext = e.value();
    fs_ascii_lower(pathext);
  }

  fs_ascii_lower(suffix);
  if(fs_trace) std::cout << "TRACE: is_exe(" << path << "):  suffix: " << suffix << " length " << suffix.length() << " pathext: " << pathext << "\n";

#if defined(__cpp_lib_string_contains) // C++23
  return pathext.contains(suffix);
#else
  return pathext.find(suffix) != std::string::npos;
#endif

#else
  const std::string cpath(path);
  return ::access(cpath.c_str(), X_OK) == 0;
#endif
}

#include <string>
#include <string_view>
#include <iostream>

#if __has_include(<ranges>)
#include <ranges>  // IWYU pragma: keep
#endif

#include <algorithm> // for std::transform, std::ranges::contains, std::replace, std::binary_search
#include <cctype> // for std::isalnum, std::toupper
#include <array>

#include "ffilesystem.h"


void fs_as_posix(std::string& path)
{
  if(fs_is_windows())
    std::replace(path.begin(), path.end(), '\\', '/');
}

std::string fs_as_posix(std::string_view path)
{
  std::string s(path);
  fs_as_posix(s);
  return s;
}


void fs_as_windows(std::string& path)
{
  if(fs_is_windows())
    std::replace(path.begin(), path.end(), '/', '\\');
}


void fs_ascii_lower(std::string& s)
{
  // for ASCII characters, not locale-aware.
  auto to_lower = [](char c) -> char {
    return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  };
#if defined(__cpp_lib_ranges)
  std::ranges::transform(s, s.begin(), to_lower);
#else
  std::transform(s.begin(), s.end(), s.begin(), to_lower);
#endif
}


void fs_ascii_upper(std::string& s)
{
  // for ASCII characters, not locale-aware.
  auto to_upper = [](char c) -> char {
      return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  };

#if defined(__cpp_lib_ranges)
  std::ranges::transform(s, s.begin(), to_upper);
#else
  std::transform(s.begin(), s.end(), s.begin(), to_upper);
#endif
}



bool
fs_is_reserved(std::string_view path)
{
  // defined as https://docs.python.org/3.13/library/os.path.html#os.path.isreserved
  if(!fs_is_windows())
    return false;

  std::string const n = fs_file_name(path);

  if(n.empty() || n == ".")
    return false;

  if(n.back() == '.' || n.back() == ' ')
    return true;

  // return true if filename contains any of :*?"<>| or ends with a space or period
  // https://docs.microsoft.com/en-us/windows/win32/fileio/naming-a-file

  if(n.find_first_of(R"(<>:"/\|?*)") != std::string::npos)
    return true;

  // don't detect ASCII control characters as reserved, since multi-byte characters may falsely trip that check

  std::string s = fs_stem(n);

  if(auto L = s.length(); L < 3 || L > 4)
    return false;

  constexpr std::array<std::string_view, 30> r = {"AUX",
    "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8", "COM9", "COM¹", "COM²", "COM³",
    "CON", "CONIN$", "CONOUT$",
    "LPT1", "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9", "LPT¹", "LPT²", "LPT³",
    "NUL", "PRN"
  };

  fs_ascii_upper(s);

#if defined(__cpp_lib_ranges_contains) // C++23
  return std::ranges::contains(r, s);
#elif defined(__cpp_lib_ranges) // C++20
  return std::ranges::binary_search(r, s);
#else
  return std::binary_search(r.begin(), r.end(), s);
#endif
}


static bool
fs_is_safe_char(const char c)
{
  // unordered_set<char>  8us
  // set<char, std::less<>>  6us
  // vector<char> 0.3us so much faster!
  // <string_view>.find_first_of is same speed as vector<char> but more readable
  constexpr std::string_view safe = "_-.~@#$%^&()[]{}+=,!";

  return std::isalnum(static_cast<unsigned char>(c)) || safe.find(c) != std::string_view::npos;
}


bool
fs_is_safe_name(std::string_view filename)
{
  // check that only shell-safe characters are in filename
  // hasn't been fully tested yet across operating systems and file systems--let us know if character(s) should be unsafe
  // does NOT check for reserved or device names
  // => filenames only, not paths
  // https://learn.microsoft.com/en-us/windows/win32/fileio/naming-a-file
  // we do not consider whitespaces, quotes, or ticks safe, as they can be confusing in shell scripts and command line usage

  if(fs_is_reserved(filename))
    return false;

#ifdef __cpp_lib_ranges // C++20
  return std::ranges::all_of(filename, fs_is_safe_char);
#else // C++11
  return std::all_of(filename.begin(), filename.end(), fs_is_safe_char);
#endif
}


bool
fs_non_ascii(std::string_view name)
{
  // check if name contains non-ASCII characters

#ifdef __cpp_lib_ranges // C++20
  return !std::ranges::all_of(name, [](int c) { return std::isprint(c); });
#else // C++11
  return !std::all_of(name.begin(), name.end(), [](int c) { return std::isprint(c); });
#endif

}


bool
fs_is_prefix(std::string_view prefix, std::string_view path)
{
  // Lexicographically, is prefix a prefix of path.

  if(prefix.empty() || path.empty())
    return false;

  std::string const pr = fs_normal(prefix);
  std::string const p = fs_normal(path);

  if (pr == p)
    return true;

  if (pr.size() > p.size())
    return false;

#ifdef __cpp_lib_starts_ends_with  // C++20
  return p.starts_with(pr);
#else
  return p.substr(0, pr.size()) == pr;
#endif

}


bool
fs_is_subdir(std::string_view subdir, std::string_view dir)
{
  // Lexicographically, is subdir a subdirectory of dir.

  if(subdir.empty() || dir.empty())
    return false;

  std::string const s = fs_normal(subdir);
  std::string const d = fs_normal(dir);

  if(s == d)
    return false;

  if(s.size() < d.size())
    return false;

#ifdef __cpp_lib_starts_ends_with  // C++20
  return s.starts_with(d);
#else
  return s.substr(0, d.size()) == d;
#endif

}

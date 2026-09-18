#include <string>
#include <string_view>

#include <system_error>
#include <cstring>

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#elif defined(_WIN32)
#include <cstdlib> // _splitpath_s, _MAX_DRIVE
#endif

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Shlwapi.h>
#endif

#include "ffilesystem.h"


char fs_filesep(){
  return fs_is_windows() ? '\\' : '/';
}


char fs_pathsep(){
  return fs_is_windows() ? ';' : ':';
}


const char* fs_devnull(){
  return fs_is_windows() ? "nul" : "/dev/null";
}


bool
fs_slash_first(std::string_view path)
{
  return !path.empty() && (path.front() == '/' || path.front() == fs_filesep());
}


std::string::size_type fs_strncpy(const char* path, char* result, const std::string::size_type buffer_size)
{
// check size before copy
  const auto L = std::strlen(path);
  if(L >= buffer_size){
    fs_print_error(path, std::make_error_code(std::errc::result_out_of_range));
    return 0;
  }

  if(L)
    #ifdef _MSC_VER
      strncpy_s(result, buffer_size, path, _TRUNCATE);
    #else
      std::strncpy(result, path, buffer_size);
    #endif

  return L;
}


bool fs_is_absolute(std::string_view path)
{
  // is path absolute?

#if defined(HAVE_CXX_FILESYSTEM) && !(defined(__MINGW32__) && !defined(__clang__) && defined(__GNUC__))
  // MinGW GCC <filesystem> .is_absolute doesn't handle UNC paths at least through GCC 15.2.0
  return Filesystem::path(path).is_absolute();
#else
  if(fs_is_windows()) {
    if(path.length() < 3)
      return false;

    // Extended-length or device path
    if(fs_win32_is_ext_path(path))
        return true;
#if defined(_WIN32)
    if(std::string cpath(path); PathIsUNCA(cpath.c_str()))
      return true;
#endif
    // Windows drive letter with slash (e.g. C: without slash is relative)
    return !(fs_root_name(path).empty()) && (path[2] == '/' || path[2] == fs_filesep());
  } else {
    return !path.empty() && path.front() == '/';
  }
#endif
}


bool fs_has_filename(std::string_view path)
{
  // does path have a filename component?
#if defined(HAVE_CXX_FILESYSTEM)
  return Filesystem::path(path).has_filename();
#else
  if (path.empty())
    return false;

  const auto i = fs_is_windows()
    ? path.find_last_of(R"(/\)")
    : path.rfind('/');

  return (i == std::string_view::npos) || (i < path.length() - 1);
#endif
}


std::string fs_file_name(std::string_view path)
{
#ifdef HAVE_CXX_FILESYSTEM
  return Filesystem::path(path).filename().string();
#else

  const auto i = fs_is_windows()
    ? path.find_last_of(R"(/\)")
    : path.rfind('/');

  return (i != std::string_view::npos)
    ? std::string(path.substr(i + 1))
    : std::string(path);
#endif
}


std::string fs_root(std::string_view path)
{
  // root_path = root_name / root_directory

#ifdef HAVE_CXX_FILESYSTEM
  return Filesystem::path(path).root_path().string();
#else
  if (std::string r = fs_root_name(path); r.empty()){
    if (fs_slash_first(path))
      return std::string(1, path.front());

    return {};
  } else if (fs_is_windows()) {
    if (path.length() >= 3 && (path[2] == '/' || path[2] == '\\')){
      return r + path[2];
    }

    return r;
  } else
    return "/";

#endif
}


std::string fs_root_name([[maybe_unused]] std::string_view path)
{

#ifdef HAVE_CXX_FILESYSTEM
  return Filesystem::path(path).root_name().string();
#elif defined(_WIN32)
  char drive[_MAX_DRIVE];
// https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/splitpath-s-wsplitpath-s
  std::string cpath(path);
  if(_splitpath_s(cpath.c_str(), drive, _MAX_DRIVE, nullptr, 0, nullptr, 0, nullptr, 0) == 0) {
    cpath = drive;
    fs_trim(cpath);
    return cpath;
  }
#endif
  return {};
}


std::string fs_stem(std::string_view path)
{
#ifdef HAVE_CXX_FILESYSTEM
  return Filesystem::path(path).filename().stem().string();
#else
  std::string r = fs_file_name(path);
  // handle special case a/..
  if (r == "..")
    return r;

  // find last dot
  const auto i = r.rfind('.');

  if (i != std::string::npos && i != 0)
    r.resize(i);

  return r;
#endif
}


std::string fs_suffix(std::string_view path)
{
#ifdef HAVE_CXX_FILESYSTEM
  return Filesystem::path(path).filename().extension().string();
#else
  const std::string p = fs_file_name(path);
  // find last dot
  const auto i = p.rfind('.');

  return (i != std::string::npos && i != 0)
    ? p.substr(i)
    : "";
#endif
}


std::string fs_join(std::string_view path, std::string_view other)
{
  // does not normalize to preserve possible symlinks
  if(path.empty() && other.empty())
    return {};
  if(path.empty())
    return std::string(other);
  if(other.empty())
    return std::string(path);

#ifdef HAVE_CXX_FILESYSTEM
  return (Filesystem::path(path) / other).generic_string();
#else
  if (other.front() == '/' || (fs_is_windows() && fs_is_absolute(other)))
    return std::string(other);

  std::string p(path);

  if(p.back() != '/')
    p.push_back('/');

  p += other;
  return p;
#endif
}


std::string fs_with_suffix(std::string_view path, std::string_view new_suffix)
{
  std::string const stem = fs_stem(path);
  // handle directory case: stem is empty
  if(stem.empty())
    return fs_join(path, new_suffix);

#ifdef HAVE_CXX_FILESYSTEM
  return Filesystem::path(path).replace_extension(new_suffix).string();
#else
  std::string const p = fs_parent(path);

  std::string r;
  if (p == ".")
    r = stem;
  else {
    r = p;
    if(p.back() != '/'){
      r.push_back('/');
    }
    r += stem;
  }

  r += new_suffix;
  return r;

#endif
}


std::string fs_lexically_normal(std::string_view path){
#ifdef HAVE_CXX_FILESYSTEM
  return Filesystem::path(path).lexically_normal().generic_string();
#else
  fs_print_error(path, std::make_error_code(std::errc::function_not_supported));
  return {};
#endif
}


std::string fs_make_preferred(std::string_view path){
#ifdef HAVE_CXX_FILESYSTEM
  return Filesystem::path(path).make_preferred().string();
#else
  fs_print_error(path, std::make_error_code(std::errc::function_not_supported));
  return {};
#endif
}

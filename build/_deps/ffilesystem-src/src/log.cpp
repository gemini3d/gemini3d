#include <system_error>

#include <iostream>

#include <string>
#include <string_view>

#if defined(_WIN32) || defined(__CYGWIN__)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <cerrno>
#endif

#if __has_include(<source_location>)
#include <source_location>
#endif

#if __has_include(<format>)
#include <format>
#endif

void fs_emit_error()
{
#if defined(_WIN32) || defined(__CYGWIN__)
  if (DWORD error = GetLastError(); error)
    std::cerr << "GetLastError: " << std::system_category().message(error) << " ";
#else
  if(errno){
    auto econd = std::generic_category().default_error_condition(errno);
    std::cerr << "errno: " << econd.message() << "\n";
  }
#endif
}


void fs_print_error(std::string_view path
#if defined(__cpp_lib_source_location)
, const std::source_location& location){
  std::string s =
#if defined(__cpp_lib_format)
    std::format("{}:{}", location.file_name(), location.line());
#else
    std::string(location.file_name()) + ":" + std::to_string(location.line());
#endif
  std::string_view fname = location.function_name();
#else
 ){
  std::string_view s;
  std::string_view fname;
#endif

  std::cerr << "ERROR: Ffilesystem: " << fname << "(" << path << ")\n" << s << "\n";

  fs_emit_error();
}


void fs_print_error(std::string_view path, const std::error_code& ec
#if defined(__cpp_lib_source_location)
, const std::source_location& location){
  std::string s =
#if defined(__cpp_lib_format)
    std::format("{}:{}", location.file_name(), location.line());
#else
    std::string(location.file_name()) + ":" + std::to_string(location.line());
#endif
  std::string_view fname = location.function_name();
#else
 ){
  std::string_view s;
  std::string_view fname;
#endif

  std::cerr << "ERROR: Ffilesystem: " << fname << "(" << path << ")\n" << s << "\n";
  if(ec)
    std::cerr << ec.message() << " " << ec.value() << "\n";

  fs_emit_error();
}


void fs_print_error(std::string_view path1, std::string_view path2
#if defined(__cpp_lib_source_location)
, const std::source_location& location){
  std::string s =
#if defined(__cpp_lib_format)
    std::format("{}:{}", location.file_name(), location.line());
#else
    std::string(location.file_name()) + ":" + std::to_string(location.line());
#endif
  std::string_view fname = location.function_name();
#else
 ){
  std::string s;
  std::string_view fname;
#endif

  std::cerr << "ERROR: Ffilesystem: " << fname << "(" << path1 <<  ", " << path2 << ")\n" << s << "\n";

  fs_emit_error();
}

void fs_print_error(std::string_view path1, std::string_view path2, const std::error_code& ec
#if defined(__cpp_lib_source_location)
, const std::source_location& location){
  std::string s =
#if defined(__cpp_lib_format)
    std::format("{}:{}", location.file_name(), location.line());
#else
    std::string(location.file_name()) + ":" + std::to_string(location.line());
#endif
  std::string_view fname = location.function_name();
#else
 ){
 std::string s;
 std::string_view fname;
#endif

  std::cerr << "ERROR: Ffilesystem: " << fname << "(" << path1 <<  ", " << path2 << ")\n" << s << "\n";

  if(ec)
    std::cerr << "C++ exception: " << ec.message() <<  " " << ec.value() << "\n";

  fs_emit_error();
}

#include <cstdlib>
#include <iostream>
#include <ostream> // for std::endl

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

#include "ffilesystem.h"


[[noreturn]]
void err(std::string_view m
#if defined(__cpp_lib_source_location)
 , const std::source_location& location = std::source_location::current()){
  std::cerr << "FAILED: " << m << " at " << location.file_name() << ":" << location.line();
#else
  ){
  std::cerr << "FAILED: " << m;
#endif

  std::cerr << " backend: " << fs_backend() << "\n";

#if defined(_WIN32) || defined(__CYGWIN__)
  DWORD error = GetLastError();
  std::cerr << std::system_category().message(error) << "\n";
#else
  if(errno){
    auto econd = std::generic_category().default_error_condition(errno);
    std::cerr << " " << econd.message();
  }
#endif

  std::cerr << std::endl;

  std::exit(EXIT_FAILURE);
}


[[noreturn]]
void skip(std::string_view m){
  std::cout << "SKIP: " << m << " backend: " << fs_backend() << "\n";
  std::exit(77);
}


void ok_msg(std::string_view m){
  std::cout << "OK: " << m << " backend: " << fs_backend() << "\n";
}

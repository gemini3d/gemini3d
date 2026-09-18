// use ffilesystem library from C++

#include <iostream>
#include <string>
#include <cstdlib>

#include "ffilesystem.h"

[[noreturn]] void err(std::string_view m){
    std::cerr << "ERROR: " << m << "\n";
    std::exit(EXIT_FAILURE);
}


int main() {

  auto cwd = fs_get_cwd();
  if (cwd.empty())
    err("current working dir not found");

  std::cout << "current working dir " << cwd << "\n";

  std::string h;
  h = fs_get_homedir();
  if (h.empty())
    err("home dir not found");

  std::cout << "home dir " << h << "\n";

  h = fs_expanduser("~");
  if (h.empty())
    err("home dir not found");

  std::cout << "expanduser('~') " << h << "\n";

  return EXIT_SUCCESS;
}

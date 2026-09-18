// Check that filesystem is capable of symbolic links with this compiler.
#include <iostream>
#include <cstdlib>
#include <string_view>

#include <filesystem>


int main(int argc, char **argv){

std::string_view exe = (argc > 1) ? argv[1] : argv[0];

auto tgt = std::filesystem::path(exe);

if(!std::filesystem::is_regular_file(tgt)){
  std::cerr << tgt << " is not a regular file\n";
  return 77;
}

auto lnk = tgt.parent_path() / "test.lnk";

if(!std::filesystem::exists(lnk)) {
  std::filesystem::create_symlink(tgt, lnk);
  std::cout << "created symlink: " << lnk << "\n";
}
auto l = std::filesystem::symlink_status(lnk);

if(!std::filesystem::exists(l)){
  std::cerr << "symlink not created: " << lnk << "\n";
  return EXIT_FAILURE;
}

if(!std::filesystem::is_symlink(l)){
  std::cerr << "is not a symlink: " << lnk << "\n";
  return EXIT_FAILURE;
}

std::cout << "OK: symlink: " << lnk << " => " << tgt << "\n";

return EXIT_SUCCESS;
}

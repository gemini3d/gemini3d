#include "ffilesystem.h"

#include <boost/ut.hpp>


int main() {
  using namespace boost::ut;

  "exe_path"_test = [] {
    std::string path = fs_exe_path();

    expect(!path.empty());
    expect(fs_exists(path));
    expect(!fs_is_dir(path));
    expect(fs_is_file(path));
    expect(fs_is_exe(path));
    expect(fs_is_readable(path));

    if(!fs_is_cygwin())
      expect(fs_is_executable_binary(path));
  };
}

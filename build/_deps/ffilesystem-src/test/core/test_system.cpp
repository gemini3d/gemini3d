#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
  using namespace boost::ut;

  "system"_test = [] {
    expect(!fs_compiler().empty()) << "unknown compiler";
    expect(neq(fs_compiler().length(), fs_get_max_path())) << "compiler has blank space";

    expect(!fs_get_username().empty());
    expect(neq(fs_get_username().length(), fs_get_max_path())) << "username has blank space";
  };
}

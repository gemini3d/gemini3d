#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

  "LibPath"_test = [] {
    std::string path = fs_lib_path();

    expect(!path.empty());
    expect(fs_exists(path));
    expect(!fs_is_dir(path));
    expect(fs_is_file(path));
};
}

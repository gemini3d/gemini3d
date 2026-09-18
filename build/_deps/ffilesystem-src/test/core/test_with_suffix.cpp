#include "ffilesystem.h"

#include <boost/ut.hpp>

using namespace boost::ut;

int main() {
  "WithSuffix"_test = [] {
    expect(fs_with_suffix("", ".h5") == ".h5");
    expect(fs_with_suffix("foo.h5", "") == "foo");
    expect(fs_with_suffix(".foo.h5", ".txt") == ".foo.txt");
    expect(fs_with_suffix(".h5", "") == ".h5");
    expect(fs_with_suffix(".h5", ".h5") == ".h5.h5");
    expect(fs_with_suffix("c:/a/hi.nc", ".h5") == "c:/a/hi.h5");
    expect(fs_with_suffix("my/file.h5", ".hdf5") == "my/file.hdf5");
    expect(fs_with_suffix("a/boo", ".h5") == "a/boo.h5");
    expect(fs_with_suffix("boo", ".h5") == "boo.h5");
    expect(fs_with_suffix("a/b/c.d/", ".hdf5") == "a/b/c.d/.hdf5");
    expect(fs_with_suffix("dir.h5/", ".hdf5") == "dir.h5/.hdf5");
    expect(fs_with_suffix("a/b/.h5", ".nc") == "a/b/.h5.nc");
};
}

#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
  using namespace boost::ut;

  "join"_test = [] {
    expect(eq(fs_join("", ""), std::string{""}));
    expect(eq(fs_join("a", ""), std::string{"a"}));
    expect(eq(fs_join("", "b"), std::string{"b"}));
    expect(eq(fs_join("a/b/", "c/"), std::string{"a/b/c/"}));
    expect(eq(fs_join("/", ""), std::string{"/"}));
    expect(eq(fs_join("", "/"), std::string{"/"}));
    expect(eq(fs_join("a", "b/"), std::string{"a/b/"}));
    expect(eq(fs_join("a/", "b/"), std::string{"a/b/"}));
    expect(eq(fs_join("a/b/../", "c/d/../"), std::string{"a/b/../c/d/../"}));
    expect(eq(fs_join("a/b", ".."), std::string{"a/b/.."}));
    expect(eq(fs_join("./a/b", "."), std::string{"./a/b/."}));
    expect(eq(fs_join("a/b", "c/d"), std::string{"a/b/c/d"}));
    expect(eq(fs_join("ab/cd", "/ef"), std::string{"/ef"}));
  };
}

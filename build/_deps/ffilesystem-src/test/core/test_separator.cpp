#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

"as_posix"_test = [] {

if(fs_is_windows()){
  std::string s = R"(a\b)";
  fs_as_posix(s);
  expect(s == "a/b");

  std::string t = R"(C:\my\path)";
  fs_as_posix(t);
  expect(t == "C:/my/path");
}

};

"path_sep"_test = [] {

if(fs_is_windows()){
  expect(eq(fs_pathsep(), ';'));
} else {
  expect(eq(fs_pathsep(), ':'));
}
};
}

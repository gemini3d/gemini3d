#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

"safe"_test = [] {

expect(!fs_is_safe_name("test/re/"));
expect(!fs_is_safe_name("test/re"));

if(fs_is_windows())
  expect(!fs_is_safe_name("hi."));
else
  expect(fs_is_safe_name("hi."));

expect(!fs_is_safe_name("hi there"));

expect(!fs_is_safe_name("日本語"));

};
}

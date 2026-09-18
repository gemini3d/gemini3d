#include <string>

#include <boost/ut.hpp>

#include "ffilesystem.h"
#include "ffilesystem_test.hpp"

int main() {
using namespace boost::ut;

"reserved_agnostic"_test = [] {
std::string ref = fs_devnull();

expect(!fs_is_reserved("."));

expect(eq(fs_normal(ref), ref));

bool b = fs_is_absolute(ref);
if (fs_is_windows())
  expect(!b);
else
  expect(b);

b = fs_is_reserved(ref);
if (fs_is_windows())
  expect(b);
else
  expect(!b);

expect(!fs_is_dir(ref));

const auto size = fs_file_size(ref);
expect(any_of{0UL, fs_unknown_size} == size);

// omitted fs_space_available() since some systems don't handle NUL /dev/null well
// e.g. Windows, macOS GCC, etc.

expect(eq(fs_expanduser(ref), ref));

expect(!fs_copy_file(ref, ref, false));

// touch is ambiguous on reserved, so omit

expect(!fs_is_symlink(ref));

};

#ifdef _WIN32
skip /
#endif
"reserved_posix"_test = [] {
  std::string ref = fs_devnull();

expect(!fs_is_exe(ref));

// NOTE: do not test
//
// create_directories(/dev/null)
// remove(/dev/null)
// create_symlink()
// set_permissionss()
//
// since if testing with "root" privileges,
// it can make the system unusable until reboot!

expect(fs_exists(ref));

expect(!fs_is_file(ref));

expect(!fs_canonical(ref, false).empty());
};
}

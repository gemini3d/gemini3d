// verify functions handle empty input OK

#include <string>
#include <string_view>

#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
  using namespace boost::ut;

  "empty_input"_test = [] {
    constexpr std::string_view e = "";

    std::string s = "";
    fs_as_posix(s);
    expect(s.empty());

    expect(fs_file_name(e).empty());
    expect(fs_stem(e).empty());
    expect(fs_join(e, e).empty());
    expect(fs_suffix(e).empty());
    expect(fs_with_suffix(e, e).empty());
    expect(!fs_is_char_device(e));
    expect(!fs_is_reserved(e));
    expect(!fs_is_symlink(e));
    expect(!fs_create_symlink(e, e));
    expect(!fs_mkdir(e));
    expect(fs_which(e).empty());
    expect(fs_root(e).empty());
    expect(!fs_exists(e));
    expect(!fs_is_absolute(e));
    expect(!fs_is_dir(e));
    expect(!fs_is_exe(e));
    expect(!fs_is_file(e));
    expect(!fs_remove(e));
    expect(fs_canonical(e, false).empty());
    expect(!fs_equivalent(e, e));
    expect(fs_expanduser(e).empty());
    expect(!fs_copy_file(e, e, false));
    expect(!fs_touch(e));

    expect(eq(fs_file_size(e), fs_unknown_size)) << "backend: " << fs_backend();

    expect(!fs_get_cwd().empty());
    expect(!fs_get_homedir().empty());

    if (!fs_is_windows()) {
      expect(eq(fs_space_available(e), fs_unknown_size)) << "backend: " << fs_backend();
      expect(!fs_set_permissions(e, 0, 0, 0));
    }
  };
}

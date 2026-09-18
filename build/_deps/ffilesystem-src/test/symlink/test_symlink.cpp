#include "ffilesystem.h"

#include <iostream>

#include <boost/ut.hpp>

namespace {

struct symlink_ctx {
  std::string cwd;
  std::string tgt;
  std::string link_file;
  std::string link_dir;
  std::string broken_link;
  std::string not_exist_tgt;
  std::string in_file;
  std::string in2;
  std::string_view nonnull_file;

  void cleanup() const {
    for (const auto& link : {tgt, link_file, link_dir, broken_link, in_file, in2}) {
      fs_remove(link);
    }
  }

  ~symlink_ctx() { cleanup(); }
};

auto setup(symlink_ctx& ctx, std::string_view name) {
  using namespace boost::ut;

  const std::string n = std::string{"TestSymlink-"} + std::string{name};
  ctx.cwd = fs_realpath(fs_get_cwd());
  // realpath is for Windows Dev Drive and Networked drives
  expect(!ctx.cwd.empty() >> fatal) << "get_cwd() should not return empty string";

  ctx.tgt = ctx.cwd + fs_filesep() + "test_" + n + "_cpp.txt";
  ctx.not_exist_tgt = ctx.cwd + fs_filesep() + "test_" + n + "_cpp.notexist";

  expect(fs_touch(ctx.tgt) >> fatal);
  expect(fs_is_file(ctx.tgt) >> fatal) << "is_file(" << ctx.tgt << ") should be true for existing regular file";

  ctx.link_file = ctx.cwd + fs_filesep() + "test_" + n + "_cpp.link";
  ctx.link_dir = ctx.cwd + fs_filesep() + "test_" + n + "_cpp.dir.link";
  ctx.broken_link = ctx.cwd + fs_filesep() + "test_" + n + "_cpp.broken";
  ctx.in_file = ctx.link_file + "-in_file";
  ctx.in2 = ctx.in_file + "-read_past_the_end_of_buffer";

  for (const auto& link : {ctx.link_file, ctx.link_dir, ctx.broken_link, ctx.in_file, ctx.in2}) {
    if (fs_is_symlink(link)) {
      std::cout << "Removing existing symlink: " << link << "\n";
      expect(fs_remove(link) >> fatal) << "Failed to remove existing symlink: " << link;
    }
  }

  expect(fs_create_symlink(ctx.tgt, ctx.link_file) >> fatal);
  std::cout << "Created symlink FILE: " << ctx.link_file << " -> " << ctx.tgt << "\n";
  expect(fs_create_symlink(ctx.cwd, ctx.link_dir) >> fatal);
  std::cout << "Created symlink DIR: " << ctx.link_dir << " -> " << ctx.cwd << "\n";

  // to create a broken link, we first create a valid link and then remove the target
  expect(fs_touch(ctx.not_exist_tgt) >> fatal);
  expect(fs_create_symlink(ctx.not_exist_tgt, ctx.broken_link) >> fatal);
  expect(fs_remove(ctx.not_exist_tgt) >> fatal);
  expect(!fs_exists(ctx.not_exist_tgt) >> fatal)
      << "exists() should be false for non-existent target: " << ctx.not_exist_tgt;
  std::cout << "Created broken symlink: " << ctx.broken_link << " -> " << ctx.not_exist_tgt << "\n";

  ctx.nonnull_file = std::string_view(ctx.in2.data(), ctx.in_file.size());
  expect(ctx.nonnull_file.back() != '\0' >> fatal);
}

} // namespace

int main() {
  using namespace boost::ut;

if (!fs_is_writable(".")) {
  skip / "create_symlink"_test = [] {};
  skip / "exists"_test = [] {};
  skip / "lexists"_test = [] {};
} else {
  "create_symlink"_test = [] {
    symlink_ctx ctx;
    setup(ctx, "CreateSymlink");

    expect(!fs_create_symlink(ctx.tgt, "")) << "create_symlink() should fail with empty link";
    expect(!fs_is_symlink(ctx.tgt) >> fatal) << "is_symlink() should be false for non-symlink file: " << ctx.tgt;
    expect(!fs_create_symlink("", ctx.link_file)) << "create_symlink() should fail with empty target";

    expect(fs_create_symlink(ctx.tgt, ctx.nonnull_file) >> fatal)
        << "create_symlink(" << ctx.tgt << ", " << ctx.nonnull_file
        << ") should succeed even if link path is not null-terminated";
    expect(fs_exists(ctx.in_file)) << "exists() should be false for non-existent file: " << ctx.in_file;

    expect(fs_is_symlink(ctx.link_file)) << "is_symlink() should be true for symlink: " << ctx.link_file;
    expect(fs_is_file(ctx.link_file))
        << "is_file(" << ctx.link_file << ") should be true for existing regular file target " << ctx.tgt;
    expect(eq(fs_read_symlink(ctx.link_file), ctx.tgt));

    // Cygwin will have /cygdrive/c and /home/ as roots
    if (!fs_is_cygwin()) {
      const std::string r = fs_canonical(ctx.link_file, true);
      expect(!r.empty() >> fatal);
      expect(eq(r.length(), ctx.tgt.length())) << r << " vs " << ctx.tgt;
      expect(fs_equivalent(r, ctx.tgt));
    }

    expect(fs_read_symlink(ctx.tgt).empty());
    expect(fs_read_symlink(ctx.not_exist_tgt).empty());
    expect(!fs_is_symlink(ctx.cwd));

    expect(fs_is_dir(ctx.link_dir)) << "is_dir(" << ctx.link_dir << ") should be true for link to existing dir";
    expect(fs_is_symlink(ctx.link_dir)) << "is_symlink() should be true for symlink: " << ctx.link_dir;
    expect(eq(fs_read_symlink(ctx.nonnull_file), ctx.tgt))
        << "read_symlink() should not read past the end of string_view buffer";
  };

  "exists"_test = [] {
    symlink_ctx ctx;
    setup(ctx, "Exists");

    expect(fs_exists(ctx.tgt)) << "exists() should be true for existing file: " << ctx.tgt;
    expect(fs_exists(ctx.link_file)) << "exists() should be true for existing symlink: " << ctx.link_file;
    expect(fs_exists(ctx.link_dir)) << "exists() should be true for existing symlink: " << ctx.link_dir;

    if (!fs_is_windows() || (fs_is_msvc() && fs_backend() == "<filesystem>")) {
      // Python os.stat() and os.lstat() handle broken symlink correctly on Windows
      expect(!fs_exists(ctx.broken_link)) << "exists() should be false for broken symlink: " << ctx.broken_link;
    }
  };

  "lexists"_test = [] {
    symlink_ctx ctx;
    setup(ctx, "Lexists");

    expect(!fs_lexists(ctx.not_exist_tgt)) << "lexists() should be false for non-existent target: " << ctx.not_exist_tgt;
    expect(!fs_lexists("")) << "lexists() should be false for empty path";
    expect(fs_lexists(ctx.tgt)) << "lexists() should be true for existing file: " << ctx.tgt;
    expect(fs_lexists(ctx.link_file)) << "lexists() should be true for existing symlink: " << ctx.link_file;
    expect(fs_lexists(ctx.link_dir)) << "lexists() should be true for existing symlink: " << ctx.link_dir;
    expect(fs_lexists(ctx.broken_link)) << "lexists() should be true for broken symlink: " << ctx.broken_link;
  };
}
}

#include "ffilesystem.h"
#include <iostream>

#include <boost/ut.hpp>

namespace {

struct symlink_ctx {
  std::string tfile;
  std::string tdir;
  std::string in_file;
  std::string tgt;
  std::string test_dir;
  std::string_view nonnull_file;

  void cleanup() const {
    for (const auto& link : {tfile, tdir}) {
      if (fs_is_symlink(link)) {
        fs_remove(link);
      }
    }
    fs_remove(tgt);
    fs_remove(test_dir);
  }

  ~symlink_ctx() { cleanup(); }
};

auto setup(symlink_ctx& ctx) {
  using namespace boost::ut;

  const std::string cwd = fs_get_cwd();
  ctx.test_dir = cwd + "/symlink_test_dir";
  expect(fs_mkdir(ctx.test_dir) >> fatal) << "Failed to create test directory: " << ctx.test_dir;

  ctx.tgt = ctx.test_dir + "/test_is_symlink_target.txt";
  expect(fs_touch(ctx.tgt) >> fatal) << "Failed to create target file: " << ctx.tgt;

  ctx.tdir = ctx.test_dir + "/link.dir";
  ctx.tfile = ctx.tdir + "/cmake_test_symlink.txt.link";

  if (fs_is_symlink(ctx.tfile)) {
    std::cout << "Removing existing test symlink file: " << ctx.tfile << "\n";
    expect(fs_remove(ctx.tfile) >> fatal) << "Failed to remove existing test symlink file: " << ctx.tfile;
  }

  if (fs_is_symlink(ctx.tdir)) {
    std::cout << "Removing existing test symlink dir: " << ctx.tdir << "\n";
    expect(fs_remove(ctx.tdir) >> fatal) << "Failed to remove existing test symlink dir: " << ctx.tdir;
  }

  expect(fs_create_symlink(ctx.test_dir, ctx.tdir) >> fatal)
      << "Failed to create symlink: " << ctx.tdir << " -> " << ctx.test_dir;
  expect(fs_is_dir(ctx.tdir) >> fatal) << ctx.tdir << " is not a directory";
  std::cout << "Created symlink DIR: " << ctx.tdir << " -> " << ctx.test_dir << "\n";

  expect(fs_create_symlink(ctx.tgt, ctx.tfile) >> fatal)
      << "Failed to create symlink: " << ctx.tfile << " -> " << ctx.tgt;
  expect(fs_is_file(ctx.tfile) >> fatal) << ctx.tfile << " is not a file after creating symlink";
  std::cout << "Created symlink FILE: " << ctx.tfile << " -> " << ctx.tgt << "\n";

  ctx.in_file = ctx.tfile + "-read_past_the_end_of_buffer";
  ctx.nonnull_file = std::string_view(ctx.in_file.data(), ctx.tfile.size());
  expect(ctx.nonnull_file.back() != '\0' >> fatal);
}

} // namespace

int main() {
  using namespace boost::ut;

  "is_symlink_file"_test = [] {
    symlink_ctx ctx;
    setup(ctx);

    expect(!fs_is_symlink("not-exist-file"));
    expect(!fs_is_symlink(""));
    expect(fs_is_symlink(ctx.tfile));
    expect(fs_is_symlink(ctx.nonnull_file)) << "is_symlink() should not read past the end of string_view buffer";
  };

  "is_symlink_dir"_test = [] {
    symlink_ctx ctx;
    setup(ctx);

    expect(fs_is_symlink(ctx.tdir));
  };
}

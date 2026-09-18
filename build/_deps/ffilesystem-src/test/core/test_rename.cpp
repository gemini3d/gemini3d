#include <string>

#include "ffilesystem.h"

#include <boost/ut.hpp>

namespace {

struct rename_ctx {
  std::string f1;
  std::string f2;
  std::string in_file;
  std::string_view nonnull_file;

  ~rename_ctx() {
    fs_remove(f2);
    fs_remove(in_file);
  }
};

void setup(rename_ctx& ctx) {
  using namespace boost::ut;

  ctx.f1 = "test_Ffs_rename.txt";
  ctx.f2 = "test_Ffs_rename2.txt";

  expect(fs_touch(ctx.f1) >> fatal);
  expect(fs_is_file(ctx.f1) >> fatal);
  if (fs_exists(ctx.f2)) {
    expect(fs_remove(ctx.f2) >> fatal);
  }

  ctx.in_file = ctx.f1 + "-invalid-memory-trailing-non-null-terminated-string_view";
  ctx.nonnull_file = std::string_view(ctx.in_file.data(), ctx.f1.size());
  expect(ctx.nonnull_file.back() != '\0' >> fatal) << "nonnull_file should not be null-terminated\n";
}

} // namespace

int main() {
  using namespace boost::ut;

if (!fs_is_writable(".")) {
  skip / "rename"_test = [] {};
} else {
  "rename"_test = [] {
    rename_ctx ctx;
    setup(ctx);

    expect(fs_rename(ctx.f1, ctx.f2));
    expect(fs_is_file(ctx.f2));

    expect(fs_rename(ctx.f2, ctx.nonnull_file) >> fatal);
    expect(!fs_is_file(ctx.in_file)) << ctx.in_file << " should have been renamed to " << ctx.f1;
    expect(fs_is_file(ctx.f1));
  };
}
}

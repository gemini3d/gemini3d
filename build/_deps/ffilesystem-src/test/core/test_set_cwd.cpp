#include "ffilesystem.h"
#include <string>

#include <boost/ut.hpp>

namespace {

struct set_cwd_ctx {
  std::string tmp;
  std::string in_dir;
  std::string_view nonnull_dir;
};

auto setup(set_cwd_ctx& ctx) {
  using namespace boost::ut;

  ctx.tmp = fs_drop_slash(fs_get_tempdir());
  expect(fs_is_dir(ctx.tmp) >> fatal);

  ctx.in_dir = "./invalid-memory-trailing-non-null-terminated-string_view";
  ctx.nonnull_dir = std::string_view(ctx.in_dir.data(), 2);
  expect(ctx.nonnull_dir.back() != '\0' >> fatal) << "nonnull_dir should not be null-terminated\n";
}

} // namespace

int main() {
  using namespace boost::ut;

  "set_cwd"_test = [] {
    set_cwd_ctx ctx;
    setup(ctx);

    expect(!fs_set_cwd(""));

    expect(fs_set_cwd(ctx.tmp) >> fatal);
    // needs to be fs_equivalent due to links, network drives, etc.
    expect(fs_equivalent(fs_get_cwd(), ctx.tmp))
        << "cwd " << fs_get_cwd() << " != " << ctx.tmp << " canonical " << fs_canonical(ctx.tmp);

    expect(fs_set_cwd(ctx.nonnull_dir) >> fatal) << "problem with non null-terminated path " << ctx.nonnull_dir;
    expect(fs_equivalent(fs_get_cwd(), ctx.tmp));
  };
}

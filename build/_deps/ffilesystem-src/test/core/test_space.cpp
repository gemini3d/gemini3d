#include "ffilesystem.h"
#include <iostream>

#include <boost/ut.hpp>

// for Windows need an invalid drive as a non-existing relative path just gives the space anyway.

namespace {

struct space_ctx {
  std::string dir;
  std::string in_dir;
  std::string_view nonnull_dir;
};

auto setup(space_ctx& ctx) {
  using namespace boost::ut;

  ctx.dir = fs_drop_slash(fs_get_tempdir());
  std::cout << "Testing space on backend " << fs_backend() << " with temp dir " << ctx.dir << "\n";

  ctx.in_dir = ctx.dir + "-invalid-memory-trailing-non-null-terminated-string_view";
  ctx.nonnull_dir = std::string_view(ctx.in_dir.data(), ctx.dir.size());
  expect(ctx.nonnull_dir.back() != '\0' >> fatal) << "nonnull_dir should not be null-terminated\n";
}

} // namespace

int main() {
  using namespace boost::ut;

  "space_available"_test = [] {
    space_ctx ctx;
    setup(ctx);

    const auto avail = fs_space_available(ctx.dir);
    expect(avail > 0 && avail < fs_unknown_size);

    expect(eq(fs_space_available("cc:/not-exist-available"), fs_unknown_size)) << "backend " << fs_backend();

    expect(neq(fs_space_available(ctx.nonnull_dir), fs_unknown_size))
        << "problem with non null-terminated path " << ctx.nonnull_dir;
  };

  "space_capacity"_test = [] {
    space_ctx ctx;
    setup(ctx);

    const auto cap = fs_space_capacity(ctx.dir);
    expect(cap > 0 && cap < fs_unknown_size);

    expect(eq(fs_space_capacity("cc:/not-exist-capacity"), fs_unknown_size)) << "backend " << fs_backend();

    expect(neq(fs_space_capacity(ctx.nonnull_dir), fs_unknown_size))
        << "problem with non null-terminated path " << ctx.nonnull_dir;
  };
}

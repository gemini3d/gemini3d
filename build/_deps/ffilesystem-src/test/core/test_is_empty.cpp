#include "ffilesystem.h"

#include <string>
#include <string_view>

#include <boost/ut.hpp>

namespace {

struct empty_ctx {
  std::string dir;
  std::string in_dir;
  std::string_view nonnull_dir;

  ~empty_ctx() {
    if (!dir.empty()) {
      fs_remove(dir);
    }
  }
};

void setup(empty_ctx& ctx) {
  using namespace boost::ut;

  ctx.dir = "ffs_is_empty_empty_dir";
  expect(fs_mkdir(ctx.dir) >> fatal);

  ctx.in_dir = ctx.dir + "/read_past_the_end_of_buffer";
  ctx.nonnull_dir = std::string_view(ctx.in_dir.data(), ctx.dir.size());
  expect(ctx.nonnull_dir.back() != '\0' >> fatal);
}

} // namespace

int main() {
  using namespace boost::ut;

  if (!fs_is_writable(".")) {
    skip / "is_empty"_test = [] {};
  } else {
  "is_empty"_test = [] {
    empty_ctx ctx;
    setup(ctx);

    expect(!fs_is_empty("."));
    expect(fs_is_empty(ctx.dir));

    expect(!fs_is_empty(ctx.dir + "/not-exist-is-empty_cpp"));

    expect(fs_is_empty(ctx.nonnull_dir))
        << "fs_is_empty() should not read past the end of string_view buffer";
  };
}
}

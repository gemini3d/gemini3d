#include "ffilesystem.h"

#include <string>

#include <boost/ut.hpp>

namespace {

struct remove_ctx {
  std::string file;
  std::string in2;
  std::string_view nonnull2;

  ~remove_ctx() {
    fs_remove(file);
  }
};

void setup(remove_ctx& ctx) {
  using namespace boost::ut;

  ctx.file = "ffs_remove_test.txt";
  expect(fs_touch(ctx.file) >> fatal);

  ctx.in2 = "./" + ctx.file;
  ctx.nonnull2 = std::string_view(ctx.in2.data(), 2);
  expect(ctx.nonnull2.back() != '\0' >> fatal) << "nonnull2 should not be null-terminated\n";
}

} // namespace

int main() {
  using namespace boost::ut;

if (!fs_is_writable(".")) {
  skip / "remove"_test = [] {};
} else {
  "remove"_test = [] {
    remove_ctx ctx;
    setup(ctx);

    expect(!fs_remove(ctx.nonnull2)) << "Failed with input not null-terminated\n";
    expect(fs_remove(ctx.file));
  };
}
}

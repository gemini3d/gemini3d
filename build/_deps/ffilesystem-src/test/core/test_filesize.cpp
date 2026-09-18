#include <string>
#include <fstream>
#include <string_view>

#include "ffilesystem.h"

#include <boost/ut.hpp>

namespace {

struct filesize_ctx {
  std::string file;
  std::string in_file;
  std::string_view nonnull_file;

  ~filesize_ctx() {
    if (!file.empty()) {
      fs_remove(file);
    }
  }
};

void setup(filesize_ctx& ctx) {
  using namespace boost::ut;

  ctx.file = "ffs_filesize_5bytes.txt";
  std::ofstream ofs(ctx.file);
  ofs << "hello";

  ctx.in_file = ctx.file + "-read_past_the_end_of_buffer";
  ctx.nonnull_file = std::string_view(ctx.in_file.data(), ctx.file.size());
  expect(ctx.nonnull_file.back() != '\0' >> fatal);
}

} // namespace

int main() {
  using namespace boost::ut;

if (!fs_is_writable(".")) {
  skip / "file_size"_test = [] {};
} else {
  "file_size"_test = [] {
    filesize_ctx ctx;
    setup(ctx);

    constexpr std::uintmax_t fsize = 5;

    expect(eq(fs_file_size(ctx.file), fsize)) << "backend " << fs_backend() << " file " << ctx.file;

    expect(eq(fs_file_size("."), fs_unknown_size)) << "backend " << fs_backend();
    expect(eq(fs_file_size("not-exist-file"), fs_unknown_size)) << "backend " << fs_backend();

    expect(eq(fs_file_size(ctx.nonnull_file), fsize)) << "backend " << fs_backend()
            << "fs_file_size() non-null-terminated path " << ctx.nonnull_file;
  };
}
}

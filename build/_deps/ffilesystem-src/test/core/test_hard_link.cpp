#include "ffilesystem.h"

#include <iostream>
#include <cstdint>
#include <string>
#include <string_view>

#include <boost/ut.hpp>

namespace {

struct device_ctx {
  std::string in2, cwd;
  std::string_view nonnull2;
};

auto setup(device_ctx& ctx) {
  using namespace boost::ut;
  ctx.cwd = fs_get_cwd();
  ctx.in2 = ctx.cwd +"/invalid-memory-trailing-non-null-terminated-string_view";
  ctx.nonnull2 = std::string_view(ctx.in2.data(), ctx.cwd.size());
  expect(ctx.nonnull2.back() != '\0' >> fatal) << "nonnull2 should not be null-terminated\n";
}

} // namespace

int main() {
  using namespace boost::ut;

  "hard_link"_test = [] {
    device_ctx ctx;
    setup(ctx);

    expect(fs_hard_link_count(".") >= 1U);

    expect(eq(fs_hard_link_count("not-exist-file"), fs_unknown_size)) << "backend " << fs_backend();

    std::cout << "the return code for errors e.g. not existing file is " << fs_unknown_size << "\n";

    expect(neq(fs_hard_link_count(ctx.nonnull2), fs_unknown_size) >> fatal)
        << "fs_hard_link_count(" << ctx.nonnull2 << ") should not return error code " << fs_unknown_size;
    expect(fs_hard_link_count(ctx.nonnull2) >= 1U);
  };

  "blk_size"_test = [] {
    device_ctx ctx;
    setup(ctx);

    auto b1 = fs_get_blksize(fs_get_tempdir());
    expect(b1 > 0U);

    b1 = fs_get_blksize(ctx.cwd);
    expect(b1 > 0U);
    expect(eq(fs_get_blksize(ctx.nonnull2), b1))
        << "fs_get_blksize(" << ctx.nonnull2 << ") should be the same as fs_get_blksize(" << fs_get_tempdir() << ")";
  };

  "device"_test = [] {
    device_ctx ctx;
    setup(ctx);

    std::string_view in = "./";
    expect(fs_st_dev(in) > static_cast<dev_t>(0));

    expect(eq(fs_st_dev(ctx.nonnull2), fs_st_dev(in)))
        << "fs_st_dev(" << ctx.nonnull2 << ") should be the same as fs_st_dev(" << in << ")";

    expect(eq(fs_st_dev("not-exist-file"), static_cast<dev_t>(0))) << "backend " << fs_backend();
  };

#if defined(_WIN32) && !defined(HAVE_GETFILEINFORMATIONBYNAME)
    skip /
#endif
  "inode"_test = [] {
    device_ctx ctx;
    setup(ctx);

    std::string_view in = "./";
    const auto inode_dot = fs_inode(in);
    expect(inode_dot > static_cast<ino_t>(0) >> fatal);

    expect(eq(fs_inode(ctx.nonnull2), inode_dot))
        << "fs_inode(" << ctx.nonnull2 << ") should be the same as fs_inode(" << in << ")";

    expect(eq(fs_inode("not-exist-file"), static_cast<ino_t>(0))) << "backend " << fs_backend();
  };
}

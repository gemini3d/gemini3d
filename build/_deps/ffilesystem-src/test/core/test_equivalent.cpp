#include <string>
#include <string_view>

#include "ffilesystem.h"

#include <boost/ut.hpp>

namespace {

struct equivalent_ctx {
  std::string self;
  std::string self_name;
  std::string in_file;
  std::string_view nonnull_file;
};

auto make_ctx(std::string_view arg0) {
  using namespace boost::ut;

    std::optional<equivalent_ctx> ctx{equivalent_ctx{}};

  ctx->self = arg0;
  ctx->self_name = fs_file_name(ctx->self);
  expect(fs_is_file(ctx->self) >> fatal) << "Test executable not found: " << ctx->self;

  ctx->in_file = ctx->self + "-read_past_the_end_of_buffer";
  ctx->nonnull_file = std::string_view(ctx->in_file.data(), ctx->self.size());
  expect(ctx->nonnull_file.back() != '\0') << "Test executable name should not end with null character";

  return ctx;
}

auto setup_inaccessible_directory(std::string& base, std::string& secret) {
  base = "ffs_equiv_inaccessible_dir";
  secret = base + "/secret";

  if (fs_exists(secret) && !fs_set_permissions(secret, 1, 1, 1)) return false;

  if (fs_exists(base) && (!fs_remove(secret) || !fs_remove(base))) return false;

  if (!fs_mkdir(secret)) return false;
  if (!fs_set_permissions(base, -1, -1, -1)) return false;

  if(fs_exists(secret)) return false;
  if(fs_is_dir(secret)) return false;
  if(fs_is_file(secret)) return false;

  return true;
}

} // namespace

int main(int argc, char** argv) {
  using namespace boost::ut;

if (!fs_equivalent(".", fs_parent(argv[0]))) {
  skip / "equivalent_filename"_test = [] {};
} else {
  "equivalent_filename"_test = [argv] {
    const auto ctx = make_ctx(argv[0]);
    expect(ctx.has_value() >> fatal);

    expect(fs_is_file(ctx->self_name) >> fatal) << "Test executable name not found in CWD: " << ctx->self_name;

    expect(fs_equivalent(ctx->self_name, "./" + ctx->self_name));
    expect(fs_equivalent(ctx->self_name, ctx->self));
    expect(fs_equivalent(ctx->self_name, ctx->nonnull_file));
    expect(fs_equivalent(ctx->self, ctx->self_name));
    expect(fs_equivalent(ctx->self, ctx->self));
  };

  "equivalent_relative"_test = [] {
    std::string s = "ffs_equiv_not-exist";
    expect(!fs_equivalent(s, s));

    std::string cwd = fs_get_cwd();
    expect(!cwd.empty() >> fatal);

    expect(fs_equivalent("..", fs_parent(cwd)));
    expect(fs_equivalent(".", "./"));
    expect(fs_equivalent(".", cwd));
    expect(!fs_equivalent("..", cwd));

    // NOTE: This can be false on networked drive or Windows Dev Drive
    // fs_equivalent(".", fs_realpath("."));
  };
}

if(fs_is_windows() || fs_is_cygwin() || fs_is_admin() || !fs_is_writable(".")) {
    skip / "equivalent_inaccessible_directory"_test = [] {};
} else {
  std::string base;
  std::string secret;

  if (setup_inaccessible_directory(base, secret)) {
    "equivalent_inaccessible_directory"_test = [base, secret] {

      expect(!fs_equivalent(secret, secret));

      expect(fs_set_permissions(base, 1, 1, 1));
      expect(fs_remove(secret));
      expect(fs_remove(base));
    };
  } else {
    skip / "equivalent_inaccessible_directory"_test = [] {};
  }
}
}

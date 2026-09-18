#include "ffilesystem.h"

#include <boost/ut.hpp>

#include <vector>
#include <string>
#include <algorithm> // for std::find

namespace {

struct which_ctx {
  std::string name;
  std::string dir;
  std::string rel;
  std::string testExe;
};

auto setup_which(which_ctx& ctx, std::string_view arg0) {
  using namespace boost::ut;

  ctx.testExe = fs_is_windows() ? "cmake.exe" : "ls";

  const std::string self{arg0};
  expect(fs_is_exe(self) >> fatal);

  ctx.name = fs_file_name(self);
  ctx.dir = fs_parent(self);
  ctx.rel = std::string("./") + ctx.name;
  expect(fs_is_file(ctx.rel) >> fatal) << ctx.rel << " does not exist";
}

auto setup_which_no_path(which_ctx& ctx) {
  using namespace boost::ut;

  ctx.testExe = fs_is_windows() ? "cmake.exe" : "ls";
  expect(fs_setenv("PATH", "") >> fatal);
}

} // namespace

int main(int argc, char** argv) {
  using namespace boost::ut;

  "which"_test = [argv] {
    which_ctx ctx;
    setup_which(ctx, argv[0]);

    expect(fs_is_absolute(fs_which(ctx.testExe)));
    expect(eq(fs_file_name(fs_which(ctx.testExe)), ctx.testExe));
    expect(fs_which(ctx.testExe, "nowhere").empty()) << ctx.testExe << " should not be found in 'nowhere'";

    expect(fs_which("/not/a/path").empty());
    expect(fs_which("").empty());

    // Build a relative path from the running program
    expect(!fs_which(ctx.rel).empty()) << ctx.rel << " is not found";

    expect(fs_which("not-exist/" + ctx.name).empty());
  };

  "which_local_dir"_test = [argv] {
    which_ctx ctx;
    setup_which(ctx, argv[0]);

    expect(fs_is_file(ctx.name) >> fatal) << ctx.name << " does not exist";
    // for Windows only: local dir is preferred
    if (fs_is_windows()) {
      expect(!fs_which(ctx.name).empty());
    } else {
      expect(fs_which(ctx.name).empty());
    }

    expect(!fs_which(ctx.name, ctx.dir).empty());
  };

  "which_no_path"_test = [] {
    which_ctx ctx;
    setup_which_no_path(ctx);

    expect(fs_which(ctx.testExe).empty());
    expect(fs_which(ctx.testExe, "nowhere").empty());
  };
}

#include <string>
#include <iostream>
#include <string_view>

#include "ffilesystem.h"

#include <boost/ut.hpp>

namespace {

struct exe_ctx {
  std::string exe;
  std::string noexe;
  std::string self;
  std::string in2;
  std::string_view nonnull2;

  ~exe_ctx() {
    if (!exe.empty()) {
      fs_remove(exe);
    }
    if (!noexe.empty()) {
      fs_remove(noexe);
    }
  }
};

void setup_ctx(exe_ctx& ctx, std::string_view test_name, std::string_view arg0) {
  using namespace boost::ut;

  const std::string n = std::string{"TestExe-"} + std::string{test_name};
  ctx.exe = "test_" + n + ".exe";
  ctx.noexe = "test_" + n + "_noexe.exe";

  ctx.self = arg0;
  if (fs_is_cygwin()) {
    ctx.self += ".exe";
  }
  expect(fs_is_file(ctx.self) >> fatal) << "Self not found: " << ctx.self;

  expect(fs_touch(ctx.exe) >> fatal);
  expect(fs_is_file(ctx.exe) >> fatal) << "Failed to create test executable file: " << ctx.exe;

  expect(fs_touch(ctx.noexe) >> fatal);
  expect(fs_is_file(ctx.noexe) >> fatal) << "Failed to create test non-executable file: " << ctx.noexe;

  expect(fs_set_permissions(ctx.exe, 0, 0, 1) >> fatal);
  expect(fs_set_permissions(ctx.noexe, 0, 0, -1) >> fatal);

  expect(!fs_get_permissions(ctx.exe).empty() >> fatal);

  ctx.in2 = ctx.self + "-invalid-memory-trailing-non-null-terminated-string_view";
  ctx.nonnull2 = std::string_view(ctx.in2.data(), ctx.self.size());
  expect(ctx.nonnull2.back() != '\0' >> fatal) << "nonnull2 should not be null-terminated\n";
}

} // namespace

int main(int argc, char** argv) {
  using namespace boost::ut;

  bool skipall = (fs_is_wsl() > 0 && fs_filesystem_type(fs_absolute(".")) == "v9fs") || !fs_is_writable(".");

  if (skipall){
    skip / "is_exe"_test = [] {};
    skip / "is_exe_bin"_test = [] {};
    skip / "perms_self"_test = [] {};
    skip / "is_exe_perms"_test = [] {};
    skip / "is_not_exe_perms"_test = [] {};
    skip / "chmod_exe"_test = [] {};
    skip / "chmod_noexe"_test = [] {};
  } else {

  "is_exe"_test = [argv] {
    exe_ctx ctx;
    setup_ctx(ctx, "IsExe", argv[0]);

    expect(!fs_is_exe(""));

    expect(!fs_is_file("not-exist"));
    expect(!fs_is_exe("not-exist"));

    expect(fs_is_exe(ctx.self)) << "Self " << ctx.self << " should be an executable file";
    expect(fs_is_exe(ctx.exe)) << "Test executable " << ctx.exe << " should be an executable file";

    expect(fs_is_exe(ctx.nonnull2)) << "problem with non null-terminated path " << ctx.nonnull2;
  };

#ifdef __CYGWIN__
  skip /
#endif
  "is_exe_bin"_test = [argv] {
    exe_ctx ctx;
    setup_ctx(ctx, "IsExeBin", argv[0]);

    // Cygwin is fussy about the full path, but it does work.
    // Cygwin wants the /cygdrive/ prefix rather than /home/username/ prefix.

    expect(fs_is_executable_binary(ctx.self)) << ctx.self << " is not executable binary";
    expect(fs_is_executable_binary(ctx.nonnull2)) << "problem with non null-terminated path " << ctx.nonnull2;

    expect(!fs_is_executable_binary(ctx.exe));
    expect(!fs_is_executable_binary(ctx.noexe));
  };

  "perms_self"_test = [argv] {
    exe_ctx ctx;
    setup_ctx(ctx, "PermsSelf", argv[0]);

    std::string p = fs_get_permissions(ctx.self);
    expect(p.size() >= 3U >> fatal);
    expect(eq(p[2], 'x'));
  };

  "is_exe_perms"_test = [argv] {
    exe_ctx ctx;
    setup_ctx(ctx, "IsExePerms", argv[0]);

    std::string p = fs_get_permissions(ctx.exe);
    expect(p.size() >= 3U >> fatal);
    expect(eq(p[2], 'x'));
  };

#if defined(_WIN32) || defined(__CYGWIN__)
  skip /
#endif
  "is_not_exe_perms"_test = [argv] {
    exe_ctx ctx;
    setup_ctx(ctx, "IsNotExePerms", argv[0]);

    expect(!fs_is_exe(ctx.noexe));

    std::string p = fs_get_permissions(ctx.noexe);
    expect(p.size() >= 3U >> fatal);
    expect(eq(p[2], '-'));
  };

  "chmod_exe"_test = [argv] {
    exe_ctx ctx;
    setup_ctx(ctx, "ChmodExe", argv[0]);

    std::string p = fs_get_permissions(ctx.exe);
    std::cout << "permissions before chmod(" << ctx.exe << ", true)  = " << p << "\n";

    expect(fs_set_permissions(ctx.exe, -1, 0, 1) >> fatal);
    // test executable even without read permission

    p = fs_get_permissions(ctx.exe);
    std::cout << "permissions after chmod(" << ctx.exe << ", true) = " << p << "\n";

    expect(fs_is_exe(ctx.exe));
    expect(p.size() >= 3U >> fatal);
    expect(eq(p[2], 'x'));
  };

#if defined(_WIN32) || defined(__CYGWIN__)
  skip /
#endif
  "chmod_noexe"_test = [argv] {
    exe_ctx ctx;
    setup_ctx(ctx, "ChmodNoExe", argv[0]);

    std::string p = fs_get_permissions(ctx.noexe);
    std::cout << "permissions before chmod(" << ctx.noexe << ", false)  = " << p << "\n";

    expect(fs_set_permissions(ctx.noexe, 0, 0, -1) >> fatal);

    p = fs_get_permissions(ctx.noexe);
    std::cout << "permissions after chmod(" << ctx.noexe << ",false) = " << p << "\n";

    expect(!fs_is_exe(ctx.noexe));
    expect(p.size() >= 3U >> fatal);
    expect(eq(p[2], '-'));
  };
}
}

#include "ffilesystem.h"

#include <string>
#include <string_view>
#include <vector>

#include <boost/ut.hpp>

namespace {

struct ondisk_ctx {
  std::string self;
  std::string self_name;
  std::string cwd;
  std::string sys_drive;
  std::string in_dir;
  std::string in_sys_dir;
  std::string_view nonnull_dir;
  std::string_view nonnull_sys_drive;
};

void setup(ondisk_ctx& ctx, std::string_view arg0) {
  using namespace boost::ut;

  ctx.cwd = fs_get_cwd();
  expect(!ctx.cwd.empty() >> fatal);

  if (fs_is_windows()) {
    auto d = fs_getenv("SystemDrive");
    expect(d.has_value() >> fatal) << "Failed to get SystemDrive";
    ctx.sys_drive = d.value();
  } else {
    ctx.sys_drive = "/";
  }

  ctx.self = std::string{arg0};
  ctx.self_name = fs_file_name(ctx.self);

  ctx.in_dir = "./invalid-memory-trailing-non-null-terminated-string_view";
  ctx.nonnull_dir = std::string_view(ctx.in_dir.data(), 2);
  expect(ctx.nonnull_dir.back() != '\0' >> fatal) << "nonnull_dir should not be null-terminated\n";

  ctx.in_sys_dir = ctx.sys_drive + "/invalid-memory-trailing-non-null-terminated-string_view";
  ctx.nonnull_sys_drive = std::string_view(ctx.sys_drive.data(), ctx.sys_drive.size());
  expect(ctx.nonnull_sys_drive.back() != '\0' >> fatal) << "nonnull_sys_drive should not be null-terminated\n";
}

} // namespace

int main(int argc, char** argv) {
  using namespace boost::ut;

  if (!fs_is_writable(".")) {
    skip / "exists"_test = [] {};
    skip / "is_dir"_test = [] {};
    skip / "is_file"_test = [] {};
    skip / "is_readable"_test = [] {};
    skip / "is_writable"_test = [] {};
    skip / "is_other"_test = [] {};
    skip / "stat_mode"_test = [] {};
    skip / "realpath"_test = [] {};
    skip / "get_mod_time"_test = [] {};
    skip / "touch"_test = [] {};
    skip / "filesystem_type"_test = [] {};
    skip / "removable"_test = [] {};
  } else {

  "exists"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    for (const auto& s : {std::string("."), std::string(".."), std::string("/"), ctx.self, ctx.self_name, ctx.cwd}) {
      expect(fs_exists(s)) << "Expected to exist: " << s;
    }

    for (const auto& s : {"ffs_exists_not-exist-file", ""}) {
      expect(!fs_exists(s)) << "Expected to not exist: " << s;
    }

    expect(fs_exists(ctx.nonnull_dir));
  };

  "is_dir"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    expect(!fs_is_dir(""));
    expect(fs_is_dir("."));
    expect(fs_is_dir(ctx.cwd));
    expect(!fs_is_dir(ctx.self));
    expect(!fs_is_dir("ffs_is_dir_not-exist-dir"));
  };

  "is_file"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    expect(fs_is_file(ctx.self));
    expect(!fs_is_file("ffs_is_file_not-exist-file"));
    expect(!fs_is_file(""));
    expect(!fs_is_file("."));
    expect(!fs_is_file(ctx.cwd));
  };

  "is_readable"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    expect(fs_is_readable("."));
    expect(fs_is_readable(ctx.self));
    expect(fs_is_readable(ctx.cwd));

    if (fs_is_windows()) {
      expect(fs_is_readable(ctx.sys_drive));

      if (fs_win32_long_paths_enabled()) {
        expect(fs_is_readable(R"(\\?\)" + ctx.sys_drive + "\\"));
      }
    }

    expect(fs_is_readable("/"));
    expect(fs_is_readable(ctx.nonnull_dir));
  };

  "is_writable"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    if (!fs_is_cygwin() && !fs_is_bsd()) {
      expect(fs_is_writable(ctx.self)) << "not every platform can detect writability of the current executable, but most can";
    }

    expect(fs_is_writable(ctx.cwd));

    if (fs_is_windows()) {
      if (fs_win32_long_paths_enabled()) {
        std::string s = R"(\\?\)" + ctx.self;
        expect(fs_is_writable(s)) << s;
      }
    } else if (!fs_is_admin() && !fs_is_cygwin()) {
      expect(!fs_is_writable("/"));
    }

    expect(fs_is_writable(ctx.nonnull_dir));
  };

  "is_other"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    expect(!fs_is_other(""));
    expect(!fs_is_other("."));
    expect(!fs_is_other(ctx.self));
    expect(!fs_is_other(ctx.cwd));
    expect(!fs_is_other("ffs_is_other_not-exist-file"));
  };

  "stat_mode"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    expect(neq(fs_st_mode(ctx.self), 0));
    expect(neq(fs_st_mode(ctx.cwd), 0));
    expect(eq(fs_st_mode("ffs_stat_mode_not-exist-file"), 0));
    expect(eq(fs_st_mode(""), 0));
    expect(neq(fs_st_mode(ctx.nonnull_dir), 0));
  };

  "realpath"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    std::string expected = fs_realpath(ctx.cwd);
    expect(!expected.empty() >> fatal);

    std::string r = fs_realpath(".");
    expect(!r.empty() >> fatal);
    expect(eq(r, expected)) << r << " vs " << expected;

    expect(fs_realpath("not-exist-realpath/b/c").empty());

    r = fs_realpath("..");
    expect(!r.empty() >> fatal);
    expect(eq(r, fs_parent(expected))) << r;

    expect(eq(fs_realpath(ctx.nonnull_dir), expected))
        << "problem with non null-terminated path " << ctx.nonnull_dir;
  };

  "get_mod_time"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    expect(fs_get_modtime(ctx.cwd) > 0);
    expect(fs_get_modtime(ctx.nonnull_dir) > 0)
        << "problem with non null-terminated path " << ctx.nonnull_dir;
  };

  "touch"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    std::string_view file = "ffs_touch_test_file";
    expect(fs_touch(file));

    auto t0 = fs_get_modtime(file);
    expect(t0 > 0);

    expect(fs_set_modtime(file));
    expect(fs_get_modtime(file) >= t0);
    expect(!fs_set_modtime("not-exist-file"));

    std::string in_file = ctx.self + "-invalid-memory-trailing-non-null-terminated-string_view";
    std::string_view nonnull_file = std::string_view(in_file.data(), ctx.self.size());
    expect(nonnull_file.back() != '\0' >> fatal) << "nonnull_file should not be null-terminated\n";

    expect(fs_touch(nonnull_file) >> fatal);
    expect(fs_is_file(nonnull_file));

    fs_remove(file);
  };

  "filesystem_type"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    std::string t = fs_filesystem_type(ctx.sys_drive);
    expect(!t.empty() >> fatal) << "Failed to get filesystem type for " << ctx.sys_drive;

    t = fs_filesystem_type(ctx.nonnull_sys_drive);
    expect(!t.empty()) << "problem with non null-terminated path " << ctx.nonnull_sys_drive;
  };

  "removable"_test = [argv] {
    ondisk_ctx ctx;
    setup(ctx, argv[0]);

    expect(!fs_is_removable(ctx.sys_drive))
        << "we assume that a CI system's system drive would not be removable";
  };
}

#ifndef _WIN32
  skip /
#endif
  "short_long"_test = [] {
    auto e = fs_getenv("PROGRAMFILES");
    expect(e.has_value() >> fatal) << "Failed to get PROGRAMFILES environment variable";
    std::string long_path = e.value();
    expect(!long_path.empty() >> fatal);

    expect(eq(fs_longname(fs_shortname(long_path)), long_path));
  };
}

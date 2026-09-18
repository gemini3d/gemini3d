#include "ffilesystem.h"

#include <string>
#include <string_view>

#include <boost/ut.hpp>
namespace {
  using namespace boost::ut;

struct absolute_ctx {
  std::string base;
  std::string ref;
  std::string cwd;
  std::string sys_drive;
};

auto make_ctx() {
  absolute_ctx ctx{};

  ctx.cwd = fs_absolute("", "");

  if (fs_is_windows()) {
    ctx.base = "j:/foo";
    ctx.ref = std::string{"j:/foo"} + fs_filesep() + "rel";
    if (const auto d = fs_getenv("SystemDrive")) {
      expect(d.has_value() >> fatal);
      ctx.sys_drive = d.value();
    }
  } else {
    ctx.base = "/foo";
    ctx.ref = "/foo/rel";
  }

  return ctx;
}

} // namespace

int main() {
  using namespace boost::ut;

  "absolute"_test = [] {
    const auto ctx = make_ctx();

    expect(!ctx.cwd.empty());
    expect(eq(fs_absolute(""), ctx.cwd));
    expect(eq(fs_absolute("", ""), ctx.cwd));
    expect(eq(fs_absolute("rel", ctx.base), ctx.ref));
    expect(eq(fs_absolute(ctx.cwd + "/rel"), ctx.cwd + "/rel"));

    // absolute("./rel") may be "/fullpath/./rel" (our method, and most <filesystem> except Windows)
    //                     or "/fullpath/rel" (Windows <filesystem>)
    // using for base "." or ".." and similar has similar ambiguity for testing.

    // relative path, empty base
    expect(eq(fs_absolute("rel", ""), ctx.cwd + fs_filesep() + "rel"));

    // empty path, relative base
    expect(eq(fs_absolute("", "rel"), ctx.cwd + fs_filesep() + "rel"));

    expect(eq(fs_absolute("日本語"), ctx.cwd + fs_filesep() + "日本語"));
    expect(eq(fs_absolute("have space"), ctx.cwd + fs_filesep() + "have space"));

    if (fs_is_windows()) {
      expect(!ctx.sys_drive.empty());
      expect(eq(fs_absolute(ctx.sys_drive + "/"), ctx.sys_drive + "/"));

      // NOTE: no, as MSYS interprets "/" totally differently depending on backend and vs MSVC
      // std::string s = fs_absolute("/");
      // fs_as_posix(s);
      // expect(eq(fs_drop_slash(s), cwd));

      if (fs_win32_long_paths_enabled()) {
        expect(eq(fs_absolute(R"(\\?\X:\anybody)"), std::string{R"(\\?\X:\anybody)"}));
        expect(eq(fs_absolute(R"(\\?\UNC\server\share)"), std::string{R"(\\?\UNC\server\share)"}));
      }
    } else {
      expect(eq(fs_absolute("/"), std::string{"/"}));
    }
  };

  "is_absolute"_test = [] {
    const auto ctx = make_ctx();

    expect(!fs_is_absolute(""));
    expect(!fs_is_absolute("日本語"));
    expect(!fs_is_absolute("some space here"));

    if (fs_is_windows()) {
      expect(!ctx.sys_drive.empty());
      expect(fs_is_absolute(ctx.sys_drive + "/"));

      expect(fs_is_absolute("J:/"));
      expect(fs_is_absolute("j:/"));
      expect(!fs_is_absolute("j:"));
      expect(!fs_is_absolute("/"));
      expect(!fs_is_absolute("/日本語"));

      if (fs_win32_long_paths_enabled()) {
        expect(fs_is_absolute(R"(\\?\)"));
        expect(fs_is_absolute(R"(\\.\)"));

        expect(fs_is_absolute(R"(\\?\C:\)"));
        expect(fs_is_absolute(R"(\\.\C:\)"));
        expect(fs_is_absolute(R"(\\?\UNC\server\share)"));
        expect(fs_is_absolute(R"(\\?\UNC\server\share\日本語)"));
        expect(fs_is_absolute(R"(\\server\share\some space here)"));
        expect(fs_is_absolute(R"(\\?\C:\some space here)"));

        const std::string unc_prefixed = R"(\\server\share\must-not-be-read)";
        const std::string_view truncated_unc(unc_prefixed.data(), 2);
        expect(eq(truncated_unc, std::string_view{R"(\\)"}) >> fatal);
        expect(!fs_is_absolute(truncated_unc));
      }
    } else {
      expect(fs_is_absolute("/"));
      expect(fs_is_absolute("/日本語"));
      expect(!fs_is_absolute("j:/"));
    }
  };
}

#include "ffilesystem.h"

#include <string>

#include <boost/ut.hpp>

namespace {

struct app_exec_ctx {
  std::string path;
};

void setup(app_exec_ctx& ctx) {
  using namespace boost::ut;

  const std::string appdir = fs_getenv("LOCALAPPDATA").value_or("") + "/Microsoft/WindowsApps";
  expect(fs_is_dir(appdir) >> fatal);

  for (const auto& exe : {"wt.exe", "winget.exe", "wsl.exe", "bash.exe"}) {
    ctx.path = fs_which(exe, appdir);
    if (!ctx.path.empty()) {
      break;
    }
  }

  expect(!ctx.path.empty() >> fatal)
    << "Failed to find app execution alias. Please make sure at least one of wt.exe, winget.exe, wsl.exe or bash.exe is installed and enabled as an app execution alias.";
}

} // namespace

int main() {
  using namespace boost::ut;

#ifndef _WIN32
  skip /
#endif
  "app_exec_alias"_test = [] {
    app_exec_ctx ctx;
    setup(ctx);

    expect(fs_is_appexec_alias(ctx.path));
  };
}

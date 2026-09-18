#include "ffilesystem.h"

#include <string_view>
#include <string>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/stat.h>  // mkfifo
#include <unistd.h> // getpid
#endif

#if __has_include(<format>)
#include <format>
#endif

#include <boost/ut.hpp>

namespace {

struct fifo_ctx {
  std::string name;
#if defined(_WIN32)
  HANDLE hPipe{INVALID_HANDLE_VALUE};
#endif

  ~fifo_ctx() {
#if defined(_WIN32)
    if (hPipe != INVALID_HANDLE_VALUE) {
      CloseHandle(hPipe);
    }
    if (!name.empty()) {
      DeleteFileA(name.c_str());
    }
#else
    if (!name.empty()) {
      unlink(name.c_str());
    }
#endif
  }
};

bool setup(fifo_ctx& ctx) {
  using namespace boost::ut;

#if defined(_WIN32)
      // https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-createnamedpipea

      // must have this path prefix or INVALID_HANDLE_VALUE results
#ifdef __cpp_lib_format  // C++20
  ctx.name = std::format(R"(\\.\pipe\test_pipe_{})", GetCurrentProcessId());
#else
  ctx.name = R"(\\.\pipe\test_pipe_)" + std::to_string(GetCurrentProcessId());
#endif
  ctx.hPipe = CreateNamedPipeA(ctx.name.data(),
                PIPE_ACCESS_DUPLEX,
                PIPE_TYPE_BYTE,
                1,
                0, 0, 0, nullptr);

  return ctx.hPipe != INVALID_HANDLE_VALUE;
#else

#ifdef __cpp_lib_format
  ctx.name = std::format("test_pipe_{}", getpid());
#else
  ctx.name = "test_pipe_" + std::to_string(getpid());
#endif

  return mkfifo(ctx.name.c_str(), 0666) == 0;
#endif
}

} // namespace

int main() {
  using namespace boost::ut;

if(fs_is_android()) {
  skip / "is_fifo"_test = [] {};
  skip / "is_file"_test = [] {};
  skip / "exists"_test = [] {};
  return 77;
}

fifo_ctx ctx;
bool ok = setup(ctx);

if (ok) {
  "is_fifo"_test = [&ctx] {
    expect(fs_is_fifo(ctx.name));
  };
} else {
  skip / "is_fifo"_test = [] {};
}

if (!ok || (fs_is_windows() && fs_backend() == "<filesystem>" &&
      (fs_is_msvc() || (fs_is_mingw() && fs_compiler().substr(0, 5) == "Clang")))) {
        skip / "is_file"_test = [] {};
} else {
  "is_file"_test = [&ctx] {
    expect(!fs_is_file(ctx.name));
  };
}

if (!ok || fs_is_windows()) {
        skip / "exists"_test = [] {};
} else {
  "exists"_test = [&ctx] {
    expect(fs_exists(ctx.name));
  };
}
}

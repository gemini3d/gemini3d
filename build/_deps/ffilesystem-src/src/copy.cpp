#if defined(__linux__)
#if !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE
#endif
#if !defined(_FILE_OFFSET_BITS)
#define _FILE_OFFSET_BITS 64
#endif
#endif

#include "ffilesystem.h"

#include <iostream> // IWYU pragma: keep
#include <system_error>

#include <string_view>
#include <algorithm>

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#else

#include <cstdlib>
#include <memory>  // for std::make_unique
#include <cerrno>  // for errno

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
// for non-Windows file loop fallback
#include <sys/types.h>  // for off_t, ssize_t
#include <sys/stat.h>
#include <unistd.h> // for read, write
#include <fcntl.h>  // for open, close

// O_CLOEXEC is not defined on all systems. We're not spawning processes,
// so it would only be meaningful possibly if we're invoked multithreaded.
#if !defined(O_CLOEXEC)
#define O_CLOEXEC  0
#endif

#endif

#if defined(FFS_DARWIN) && __has_include(<copyfile.h>)
// macOS 10.12 or later
#define HAVE_MACOS_COPYFILE
#include <copyfile.h>
#endif

#endif  // HAVE_CXX_FILESYSTEM


#if !defined(HAVE_CXX_FILESYSTEM) && !defined(_WIN32)

// for RAII file descriptor management
struct fd_handle {
  int fd{-1};
  explicit fd_handle(int v = -1) : fd(v) {}
  ~fd_handle() { reset(); }

  fd_handle(const fd_handle&) = delete;
  fd_handle& operator=(const fd_handle&) = delete;
  fd_handle(fd_handle&& o) noexcept : fd(o.fd) { o.fd = -1; }
  fd_handle& operator=(fd_handle&& o) noexcept {
    if (this != &o) { reset(); fd = o.fd; o.fd = -1; }
    return *this;
  }

  int close() {
    int rv = 0;
    if (fd != -1) {
      rv = ::close(fd);
      fd = -1;
    }
    return rv;
  }

  int get() const { return fd; }
  explicit operator bool() const { return fd != -1; }

  int release() {
    int tmp = fd; fd = -1;
    return tmp;
  }

  void reset(int v = -1) {
    if (fd != -1) ::close(fd);
    fd = v;
  }
};


static off_t fs_copy_loop(int const rid, int const wid, off_t const len)
{
  // copy a file in chunks
  off_t r = len;

  if (fs_trace) std::cout << "TRACE::ffilesystem:copy_file: using plain file buffer read / write\n";

  constexpr size_t bufferSize = 16384;
  auto buf = std::make_unique<char[]>(bufferSize);

  while (r > 0) {
    const size_t to_read = static_cast<size_t>(std::min<off_t>(r, static_cast<off_t>(bufferSize)));

    ssize_t ret;
    do {
      ret = ::read(rid, buf.get(), to_read);
    } while (ret < 0 && errno == EINTR);

    if (ret <= 0)
      break;

    ssize_t written = 0;
    while (written < ret) {
      ssize_t w;
      do {
        w = ::write(wid, buf.get() + written, static_cast<size_t>(ret - written));
      } while (w < 0 && errno == EINTR);

      if (w <= 0) {
        ret = -1;
        break;
      }
      written += w;
    }

    if (ret < 0)
      break;

    r -= written;
  }

  return r;
}


bool fs_copy_file_range_or_loop(std::string_view source, std::string_view dest, bool overwrite)
{
  // copy a file in chunks
  const std::string src(source);
  fd_handle rid(::open(src.c_str(), O_RDONLY | O_CLOEXEC));
  if (!rid)
    return false;

  // leave fstat here to avoid source file race condition
  struct stat s;
  if (::fstat(rid.get(), &s) == -1)
    return false;

  const off_t len = s.st_size;

  auto opt = O_CREAT | O_WRONLY | O_TRUNC | O_CLOEXEC;
  if(!overwrite)
    opt |= O_EXCL;

// https://linux.die.net/man/3/open
  const std::string dst(dest);
  fd_handle wid(::open(dst.c_str(), opt, s.st_mode));
  if (!wid)
    return false;

  off_t r = len;
  off_t ret = 0;
  int cfr_errno = 0;

#if defined(ffilesystem_HAVE_COPY_FILE_RANGE)
    // https://man.freebsd.org/cgi/man.cgi?copy_file_range(2)
    // https://man7.org/linux/man-pages/man2/copy_file_range.2.html

  std::string const fst = fs_filesystem_type(source);

  bool useloop = fst == "debugfs" ||
                 fst == "procfs" ||
                 fst == "sysfs" ||
                 fst == "tracefs";

  if (!useloop) {
    if (fs_trace) std::cout << "TRACE::ffilesystem:copy_file: using copy_file_range\n";

    while (r > 0) {
      ssize_t cfr;
      do {
        cfr = ::copy_file_range(rid.get(), nullptr, wid.get(), nullptr, r, 0);
      } while (cfr < 0 && errno == EINTR);

      ret = cfr;
      if (cfr <= 0) {
        if (cfr < 0) cfr_errno = errno;
        break;
      }

      r -= cfr;
    }
  }
#endif

  // https://github.com/boostorg/filesystem/issues/184
  if (r != 0 || (ret < 0 && (cfr_errno == EINVAL || cfr_errno == EOPNOTSUPP))) {
    if (fs_trace) std::cout << "TRACE::ffilesystem:copy_file: falling back to fs_copy_loop (r=" << r << ", errno=" << (cfr_errno ? cfr_errno : errno) << ")\n";
    r = fs_copy_loop(rid.get(), wid.get(), r);
  }

  int wc = wid.close();
  int rc = rid.close();
  return r == 0 && wc == 0 && rc == 0;
}
#endif


bool fs_copy_file(std::string_view source, std::string_view dest, bool overwrite)
{
  // copy a single file -- symlinks are followed

  std::error_code ec;

#ifdef HAVE_CXX_FILESYSTEM
  auto opt = Filesystem::copy_options::none;
  if (overwrite)
    opt = Filesystem::copy_options::overwrite_existing;
// WORKAROUND: Windows MinGW GCC 11..13, Intel oneAPI Linux: bug with overwrite_existing failing on overwrite

  if(overwrite && fs_is_file(dest) && !fs_remove(dest))
    fs_print_error(dest, std::make_error_code(std::errc::io_error));

  if(Filesystem::copy_file(source, dest, opt, ec) && !ec)
    return true;

#elif defined(_WIN32)
  DWORD opts = 0;
  if(!overwrite)
    opts |= COPY_FILE_FAIL_IF_EXISTS;

  // https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-copyfileexw
  // preserves source file attributes

  if(CopyFileExW(fs_win32_to_wide(source).c_str(),
                 fs_win32_to_wide(dest).c_str(),
                 nullptr, nullptr, FALSE, opts) != 0)
    return true;

#elif defined(HAVE_MACOS_COPYFILE)
  // https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/copyfile.3.html
  // preserves source file attributes
  auto opt = COPYFILE_ALL;
  if (!overwrite)
    opt |= COPYFILE_EXCL;

  const std::string src = std::string(source);
  const std::string dst = std::string(dest);
  if(::copyfile(src.c_str(), dst.c_str(), nullptr, opt) == 0)
    return true;
#else

  if(fs_copy_file_range_or_loop(source, dest, overwrite))
    return true;

#endif

  fs_print_error(source, dest, ec);
  return false;

}

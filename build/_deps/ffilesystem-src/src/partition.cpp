#include <string>
#include <string_view>

#include <iostream>
#include <system_error>

#include "ffilesystem.h"

#if defined(__linux__) && (__has_include(<linux/magic.h>) || __has_include("linux/magic.h"))
// GCC < 10 doesn't detect <linux/magic.h>
// IWYU pragma: no_include <sys/statfs.h>
#define HAVE_LINUX_MAGIC_H
#include <sys/vfs.h> // IWYU pragma: keep
#include <linux/magic.h>
// https://github.com/torvalds/linux/blob/master/include/uapi/linux/magic.h
#elif defined(FFS_DARWIN) || defined(FFS_BSD)
#include <sys/param.h>
#include <sys/mount.h>
#elif defined(_WIN32) || defined(__CYGWIN__)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif



#ifdef HAVE_LINUX_MAGIC_H
static inline std::string fs_type_linux(std::string_view path)
{
  struct statfs s;

  const std::string cpath(path);

  if(statfs(cpath.c_str(), &s)) {
    fs_print_error(path);
    return {};
  }

#ifndef BTRFS_SUPER_MAGIC
#define BTRFS_SUPER_MAGIC 0x9123683E
#endif
#ifndef DEBUGFS_MAGIC
#define DEBUGFS_MAGIC 0x64626720
#endif
#ifndef FUSE_SUPER_MAGIC
#define FUSE_SUPER_MAGIC 0x65735546
#endif
#ifndef EXFAT_SUPER_MAGIC
#define EXFAT_SUPER_MAGIC 0x2011BAB0
#endif
#ifndef EROFS_SUPER_MAGIC_V1
#define EROFS_SUPER_MAGIC_V1 0xE0F5E1E2
#endif
#ifndef F2FS_SUPER_MAGIC
#define F2FS_SUPER_MAGIC 0xF2F52010
#endif
#ifndef PROC_SUPER_MAGIC
#define PROC_SUPER_MAGIC 0x9FA0
#endif
#ifndef SYSFS_MAGIC
#define SYSFS_MAGIC 0x62656572
#endif
#ifndef TRACEFS_MAGIC
#define TRACEFS_MAGIC 0x74726163
#endif
#ifndef UDF_SUPER_MAGIC
#define UDF_SUPER_MAGIC 0x15013346
#endif
#ifndef XFS_SUPER_MAGIC
#define XFS_SUPER_MAGIC 0x58465342
#endif

  switch (s.f_type) {
    case BTRFS_SUPER_MAGIC: return "btrfs";
    case DEBUGFS_MAGIC: return "debugfs";
    case EXT4_SUPER_MAGIC: return "ext4";
    case EXFAT_SUPER_MAGIC: return "exfat";
    case EROFS_SUPER_MAGIC_V1: return "erofs";
    case F2FS_SUPER_MAGIC: return "f2fs";
    case FUSE_SUPER_MAGIC: return "fuse";
    case NFS_SUPER_MAGIC: return "nfs";
    case PROC_SUPER_MAGIC: return "procfs";
    case SQUASHFS_MAGIC: return "squashfs";
    case SYSFS_MAGIC: return "sysfs";
    case TMPFS_MAGIC: return "tmpfs";
    case TRACEFS_MAGIC: return "tracefs";
    case UDF_SUPER_MAGIC: return "udf";
    case V9FS_MAGIC: return "v9fs";
    // used for WSL
    // https://devblogs.microsoft.com/commandline/whats-new-for-wsl-in-windows-10-version-1903/
    case XFS_SUPER_MAGIC: return "xfs";

    default:
      std::cerr << "ERROR:fs_filesystem_type " << path << " unknown type ID: " << s.f_type << "\n";
      return {};
  }
}
#endif


std::string fs_filesystem_type(std::string_view path)
{
  // return name of filesystem type if known

  std::error_code ec;

#if defined(_WIN32) || defined(__CYGWIN__)

  std::string r(path);

  // Cygwin: assume user input Cygwin path root directly.
  if(!fs_is_cygwin()){

    r = fs_root_name(r);
    if(r.empty())
      return {};

    // GetVolumeInformationA requires a trailing backslash
    r.push_back('\\');
  }

  if(fs_trace) std::cout << "TRACE:filesystem_type(" << path << ") root: " << r << "\n";

  // https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getvolumeinformationa
  std::string name;
  name.resize(MAX_PATH+1);

  if(GetVolumeInformationA(r.c_str(), nullptr, 0, nullptr, nullptr, nullptr, name.data(), static_cast<DWORD>(name.size()))) {
    fs_trim(name);
    return name;
  }

#elif defined(__linux__)
# ifdef HAVE_LINUX_MAGIC_H
  return fs_type_linux(path);
# else
  ec = std::make_error_code(std::errc::function_not_supported);
# endif
#elif defined(FFS_DARWIN) || defined(FFS_BSD)
  struct statfs s;
  const std::string cpath(path);

  if(!::statfs(cpath.c_str(), &s))
    return s.f_fstypename;
#else
  ec = std::make_error_code(std::errc::function_not_supported);
#endif

  fs_print_error(path, ec);
  return {};
}

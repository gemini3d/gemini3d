#include <string>

#include "ffilesystem.h"

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#else

#include <sys/stat.h>
#if __has_include(<format>)
#include <format> // IWYU pragma: keep
#endif

#endif

#if defined(__linux__)
#include <fstream>  // for std::ifstream
#include <sys/sysmacros.h> // for minor, major
#endif

#if defined(FFS_DARWIN)
#include <CoreFoundation/CoreFoundation.h>
#include <DiskArbitration/DiskArbitration.h>
#endif


bool
fs_is_removable(std::string_view path)
{
  // is path a removable device like a USB stick or SD card or CD-ROM, DVD, Blu-ray
  // not a fixed disk like a hard drive or SSD
#if defined(_WIN32)
  // https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-getdrivetypew

  UINT t = GetDriveTypeW(fs_win32_to_wide(fs_root(path)).c_str());
  switch (t)
  {
    case DRIVE_REMOVABLE: case DRIVE_CDROM:
      return true;
    case DRIVE_FIXED: case DRIVE_REMOTE: case DRIVE_RAMDISK:
      return false;
    case DRIVE_UNKNOWN: case DRIVE_NO_ROOT_DIR:
      fs_print_error(path, std::make_error_code(std::errc::no_such_device));
      return false;
    default:
      return false;
  }
#elif defined(__linux__)
  // Linux: find the device and check /sys/block/*/removable == 1
  // this is for local devices only, not network mounts
  // for more general/robust solution, consider libudev.

  std::string dev;
  const std::string cpath(path);

  // https://man7.org/linux/man-pages/man2/stat.2.html

  if (struct stat s; ::stat(cpath.c_str(), &s) == 0) {

  // don't call as ::major or ::minor because some platforms e.g. Android have them as a macro
#if defined(__cpp_lib_format)  // C++20
    dev = std::format("/sys/dev/block/{}:{}/removable", major(s.st_dev), minor(s.st_dev));
#else
    dev = "/sys/dev/block/" + std::to_string(major(s.st_dev)) + ":" + std::to_string(minor(s.st_dev)) + "/removable";
#endif
  } else {
    fs_print_error(path);
    return false;
  }

  {
    std::ifstream ifs(dev);
    if (ifs) {
      char c;
      if (ifs.get(c))
        return c == '1';
    } else if (errno == EACCES || errno == ENOENT) {
      // Android / SELinux may deny access to sysfs; virtual / network filesystems
      // won't have a sysfs entry. In both cases assume not removable.
      return false;
    }
  }

  fs_print_error(dev);
  return false;

#elif defined(FFS_DARWIN)
  // macOS DiskArbitration framework
  // requires linking with -framework DiskArbitration -framework CoreFoundation

  // https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man2/stat.2.html

  const std::string cpath(path);

  struct stat s;
  if (::stat(cpath.c_str(), &s) != 0) {
    fs_print_error(path);
    return false;
  }

   // Construct BSD device name (e.g., "disk0s1")
  std::string bsdName;

  // don't call as ::major or ::minor because some platforms e.g. Android have them as a macro

#if defined(__cpp_lib_format)  // C++20
  bsdName = std::format("disk{}s{}", major(s.st_dev), minor(s.st_dev));
#else
  bsdName = "disk" + std::to_string(major(s.st_dev)) + "s" + std::to_string(minor(s.st_dev));
#endif

  DASessionRef session = DASessionCreate(kCFAllocatorDefault);
  if (!session) {
    fs_print_error(path);
    return false;
  }

  bool ok = false;
  if (DADiskRef disk = DADiskCreateFromBSDName(kCFAllocatorDefault, session, bsdName.c_str()); disk) {
    if (CFDictionaryRef desc = DADiskCopyDescription(disk); desc) {
      if (CFBooleanRef r = static_cast<CFBooleanRef>(CFDictionaryGetValue(desc, kDADiskDescriptionMediaRemovableKey)); r) {
        ok = CFBooleanGetValue(r);
      }
      CFRelease(desc);
    }
    CFRelease(disk);
  }

  CFRelease(session);

  return ok;

#else
  fs_print_error(path, std::make_error_code(std::errc::function_not_supported));
  return false;
#endif

}

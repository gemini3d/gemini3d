#include "ffilesystem.h"

#if defined(_WIN32) || defined(__CYGWIN__)
# define WIN32_LEAN_AND_MEAN
# ifndef NOMINMAX
#  define NOMINMAX
# endif
# include <windows.h>
#else
# include <unistd.h> // for sysconf

# if defined(FFS_DARWIN)
#  include <mach/mach.h>
#  include <sys/sysctl.h>
# elif defined(__linux__)
#  include <sys/sysinfo.h>
# elif defined(FFS_BSD)
#  include <sys/types.h>
#  include <sys/sysctl.h>
# endif
#endif

#include <limits>



unsigned long long fs_total_sys_memory()
{
// https://stackoverflow.com/a/2513561

#if defined(_WIN32) || defined(__CYGWIN__)

  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if(GlobalMemoryStatusEx(&status))
    return status.ullTotalPhys;

  fs_print_error("", "GlobalMemoryStatusEx");

#else

  long pages = ::sysconf(_SC_PHYS_PAGES);

  if (pages > 0) {
    long page_size = ::sysconf(_SC_PAGE_SIZE);

    if (page_size > 0)
      return static_cast<unsigned long long>(pages) * page_size;
  }

  fs_print_error("", "sysconf");

#endif

  return 0;
}


unsigned long long fs_get_free_memory()
{
  // https://github.com/ninja-build/ninja/pull/2605

#if defined(_WIN32) || defined(__CYGWIN__)
  thread_local static unsigned long long committed_idle = std::numeric_limits<unsigned long long>::max();
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (!GlobalMemoryStatusEx(&status)) {
    fs_print_error("", "GlobalMemoryStatusEx");
    return std::numeric_limits<unsigned long long>::max();
  }
  const unsigned long long committed = (status.ullTotalPageFile - status.ullAvailPageFile);
  // since system use committed memory normally, store the smallest amount we have seen to guess how much
  // paging is not related to our program's needs
  committed_idle = (committed_idle < committed) ? committed_idle : committed;

  // this checks for wraparound
  return (status.ullAvailPhys > (committed - committed_idle))
    ? status.ullAvailPhys - (committed - committed_idle)
    : std::numeric_limits<unsigned long long>::max();

#elif defined(FFS_DARWIN)
  thread_local static unsigned long long swapped_idle = std::numeric_limits<unsigned long long>::max();
  vm_size_t page_size;
  vm_statistics64_data_t vm_stats;

  mach_port_t port = ::mach_host_self();
  mach_msg_type_number_t count = sizeof(vm_stats) / sizeof(natural_t);

  size_t swap_stats_size = sizeof(struct xsw_usage);
  struct xsw_usage swap_stats;

  int ctl[] = {CTL_VM, VM_SWAPUSAGE};

  int result = ::host_page_size(port, &page_size);
  result |= ::host_statistics64(port, HOST_VM_INFO, reinterpret_cast<host_info64_t>(&vm_stats), &count);
  result |= ::sysctl(ctl, 2, &swap_stats, &swap_stats_size, nullptr, 0);

  if (KERN_SUCCESS != result)
    return std::numeric_limits<unsigned long long>::max();

  // information not available

  // inactive memory that is marked to be moved to swap or is fs cache and should be considered as free
  unsigned long long free_memory = (vm_stats.free_count + vm_stats.inactive_count ) * page_size;
  // since inactive memory can be moved to the swap this value will be inexact.
  unsigned long long used_swap = swap_stats.xsu_used;
  swapped_idle = (used_swap < swapped_idle) ? used_swap : swapped_idle;
  return free_memory - (used_swap < swapped_idle ? 0 : (used_swap - swapped_idle));
#elif defined(__linux__)
  thread_local static unsigned long long swapped_idle = std::numeric_limits<unsigned long long>::max();
  struct sysinfo infos;

  if(::sysinfo(&infos) == 0 ) {
    const unsigned long long swapped = (infos.totalswap - infos.freeswap);
    // since system use committed memory normally, store the smallest amount we have seen
    swapped_idle = (swapped_idle < swapped) ? swapped_idle : swapped;
    return infos.freeram - (swapped < swapped_idle ? 0 : (swapped - swapped_idle));
  }
#elif defined(FFS_BSD)
  // free + inactive + cache (like Linux MemAvailable / "usable")
  unsigned int free_count = 0, inactive_count = 0, cache_count = 0;
  std::size_t len = sizeof(free_count);
  long pagesize = sysconf(_SC_PAGESIZE);

  // https://man.freebsd.org/cgi/man.cgi?query=sysctlbyname
  if (sysctlbyname("vm.stats.vm.v_free_count", &free_count, &len, nullptr, 0) == 0) {
    len = sizeof(inactive_count);
    sysctlbyname("vm.stats.vm.v_inactive_count", &inactive_count, &len, nullptr, 0);

    len = sizeof(cache_count);
    sysctlbyname("vm.stats.vm.v_cache_count", &cache_count, &len, nullptr, 0);

    return static_cast<unsigned long long>(free_count + inactive_count + cache_count) * pagesize;
  }
#endif

  return std::numeric_limits<unsigned long long>::max();

}

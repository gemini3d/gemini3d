#include <windows.h>

int cpu_count_c(void){

  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  int nCPU = sysinfo.dwNumberOfProcessors;

  if (nCPU >= 2) nCPU /= 2;

  return nCPU;
}

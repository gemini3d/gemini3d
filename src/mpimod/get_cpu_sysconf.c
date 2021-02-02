#include <unistd.h>


int cpu_count_c(void){
  int nCPU = sysconf(_SC_NPROCESSORS_ONLN);

  if (nCPU >= 2) nCPU /= 2;  // assume hyperthreading

  return nCPU;
}

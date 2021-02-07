#include <stdio.h>
#include <stdlib.h>
#include "get_cpu_hwloc.h"

int main(int argc, char *argv[]) {

  if(argc < 2){
    fprintf(stderr, "please input expected CPU physical core count\n");
    return 1;
  }

  int N = atoi(argv[1]);

  int Ncore = cpu_count_c();

  if (Ncore < 1){
    fprintf(stderr, "hwloc did not detect CPU count\n");
    return 1;
  }

  if (Ncore < 2){
    fprintf(stderr, "hwloc may not have detected CPU count\n");
    return 1;
  }

  if (Ncore == N/2){
    fprintf(stderr, "did you input logical CPU count %d?  hwloc reports: %d\n", N, Ncore);
    return 2;
  }

  if (Ncore != N){
    fprintf(stderr, "CPU count mismatch:  hwloc: %d  expected: %d", Ncore, N);
    return 2;
  }

  printf("OK: hwloc CPU count: %d", Ncore);

  return 0;

}

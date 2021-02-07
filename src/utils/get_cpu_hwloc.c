// inspired by: https://stackoverflow.com/a/29414957

#include <hwloc.h>
#include <stdio.h>

int cpu_count_c(void){

  hwloc_topology_t sTopology;

  if (hwloc_topology_init(&sTopology) != 0){
    fprintf(stderr, "hwloc: could not init topology\n");
    return -1;
  }
  if (hwloc_topology_load(sTopology) != 0){
    fprintf(stderr, "hwloc: could not load topology\n");
    return -1;
  }
// https://www.open-mpi.org/projects/hwloc/doc/v2.4.0/a00154.php#gacd37bb612667dc437d66bfb175a8dc55
  int nCore = hwloc_get_nbobjs_by_type(sTopology, HWLOC_OBJ_CORE);
  if (nCore < 1) {
    // assume hyperthreading / 2
    nCore = hwloc_get_nbobjs_by_type(sTopology, HWLOC_OBJ_PU) / 2;
    printf("hwloc: fallback to PU count/2: %d CORE count not available\n", nCore);
  }

  hwloc_topology_destroy(sTopology);

  return nCore;

}

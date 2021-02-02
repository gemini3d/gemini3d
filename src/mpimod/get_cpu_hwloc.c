// https://stackoverflow.com/a/29414957

#include <hwloc.h>

int cpu_count_c(void){

  int nPhysicalProcessorCount = 0;
  hwloc_topology_t sTopology;

  if (hwloc_topology_init(&sTopology) == 0 && hwloc_topology_load(sTopology) == 0)
  {
      nPhysicalProcessorCount = hwloc_get_nbobjs_by_type(sTopology, HWLOC_OBJ_CORE);

      hwloc_topology_destroy(sTopology);
  }

  return nPhysicalProcessorCount;

}

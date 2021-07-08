# --- test parallelism
include(ProcessorCount)

function(cmake_cpu_count)
# on ARM e.g. Raspberry Pi, the usually reliable cmake_host_system_info gives 1 instead of true count
# fallback to less reliable ProcessorCount which does work on Raspberry Pi.

ProcessorCount(_ncount)
cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)

if(Ncpu EQUAL 1 AND _ncount GREATER 0)
  set(Ncpu ${_ncount})
endif()

set(Ncpu ${Ncpu} PARENT_SCOPE)

endfunction(cmake_cpu_count)

if(DEFINED ENV{CTEST_PARALLEL_LEVEL})
  set(Ncpu $ENV{CTEST_PARALLEL_LEVEL})
else()
  cmake_cpu_count()
endif()
message(STATUS "using Ncpu = ${Ncpu}")

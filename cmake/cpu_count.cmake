# --- test parallelism
function(cmake_cpu_count)
# on ARM e.g. Raspberry Pi, the usually reliable cmake_host_system_info gives 1 instead of true count

cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)

if(Ncpu LESS 1)
  set(Ncpu 1)
endif()

set(Ncpu ${Ncpu} PARENT_SCOPE)

endfunction(cmake_cpu_count)

if(DEFINED ENV{CTEST_PARALLEL_LEVEL})
  set(Ncpu $ENV{CTEST_PARALLEL_LEVEL})
else()
  cmake_cpu_count()
endif()
message(STATUS "parallel CTest Ncpu = ${Ncpu}")

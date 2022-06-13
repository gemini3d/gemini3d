# Ref: https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/ProcessorCount.cmake
# https://gitlab.kitware.com/search?group_id=415&project_id=541&repository_ref=master&scope=blobs&search=GetNumberOfPhysicalCPU&snippets=false
# without hwloc, hyperthreading may give incorrect result

add_library(hwloc_ifc OBJECT hwloc_ifc.f90)
target_link_libraries(hwloc_ifc PRIVATE hwloc_c)

# separate C lib to avoid compiler warnings from C <-> Fortran add_compile_options etc.
add_library(hwloc_c OBJECT)

if(NOT mpi)
  target_sources(hwloc_c PRIVATE dummy_cpu.c)
  return()
endif()


if(HWLOC_FOUND)
  target_sources(hwloc_c PRIVATE get_cpu_hwloc.c)
  target_compile_definitions(hwloc_c PRIVATE $<$<BOOL:${MSVC}>:_CRT_SECURE_NO_WARNINGS>)
  target_include_directories(hwloc_c PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )
  target_link_libraries(hwloc_c PRIVATE HWLOC::HWLOC)

  return()
endif()

check_include_file(unistd.h HAVE_UNISTD_H)
if(HAVE_UNISTD_H)
  check_symbol_exists(_SC_NPROCESSORS_ONLN unistd.h HAVE_GNU_COUNT)
endif()

if(HAVE_GNU_COUNT)
  # https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/sysconf.3.html
  target_sources(hwloc_c PRIVATE get_cpu_sysconf.c)
elseif(WIN32)
  target_sources(hwloc_c PRIVATE get_cpu_win.c)
else()
  # this should be almost never encountered
  message(STATUS "did not find libhwloc or sysconf. gemini3d.run -n Ncpu option will be needed for those using it.")
  target_sources(hwloc_c PRIVATE dummy_cpu.c)
endif()

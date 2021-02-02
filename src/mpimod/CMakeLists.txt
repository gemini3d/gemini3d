add_library(autogrid autogrid.f90)

if(realbits EQUAL 32)
  set(mpi_realprec mpi_real)
else()
  set(mpi_realprec mpi_double_precision)
endif()

configure_file(mpimod.in.f90 mpimod.f90 @ONLY)

add_library(mpimod
  ${CMAKE_CURRENT_BINARY_DIR}/mpimod.f90
  mpisend.f90 mpirecv.f90 mpihalo.f90)
# not OBJECT to avoid CMake 3.13 trouble
target_link_libraries(mpimod PRIVATE autogrid const MPI::MPI_Fortran)
if(NOT mpi)
  # the reason for this is the MPI-1 Fortran 77 linking in MUMPS
  # this seemed to be an easy way to avoid symbol clashes
  target_link_libraries(mpimod PRIVATE MUMPS::mpiseq)
endif()

# --- Fortran front end
add_executable(gemini3d.run runner.f90 gemini_runner.f90 get_cpu.f90)
target_link_libraries(gemini3d.run PRIVATE const config pathlib timeutils autogrid reader h5fortran::h5fortran nc4fortran::nc4fortran)
set_target_properties(gemini3d.run PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set_target_properties(gemini3d.run PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(gemini3d.run INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

# Ref: https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/ProcessorCount.cmake
# https://gitlab.kitware.com/search?group_id=415&project_id=541&repository_ref=master&scope=blobs&search=GetNumberOfPhysicalCPU&snippets=false
# shortcomings: hyperthreading may give incorrect result
if(mpi)
  if(WIN32)
    # Windows
    target_sources(gemini3d.run PRIVATE get_cpu_win.c)
  else()
    # MacOS and Linux
    find_path(hwloc_inc NAMES hwloc.h)
    find_library(hwloc NAMES hwloc)
    if(hwloc_inc AND hwloc)
      target_sources(gemini3d.run PRIVATE get_cpu_hwloc.c)
      target_include_directories(gemini3d.run PRIVATE ${hwloc_inc})
      target_link_libraries(gemini3d.run PRIVATE ${hwloc})
    else()
      target_sources(gemini3d.run PRIVATE get_cpu_sysconf.c)
    endif()
  endif()
else(mpi)
  target_sources(gemini3d.run PRIVATE dummy_cpu.c)
endif(mpi)

#-----

if(${PROJECT_NAME}_BUILD_TESTING)

add_executable(test_excess_mpi test_excess_mpi.f90)
target_link_libraries(test_excess_mpi PRIVATE autogrid)
add_test(NAME unit:mpi:excessCPU COMMAND $<TARGET_FILE:test_excess_mpi>)
set_tests_properties(unit:mpi:excessCPU PROPERTIES TIMEOUT 10)

endif()

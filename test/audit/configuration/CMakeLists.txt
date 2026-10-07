cmake_minimum_required(VERSION 3.26)
project(test_configuration LANGUAGES NONE)
enable_testing()
include(FetchContent)

get_filename_component(gemini_source "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
include("${gemini_source}/test/mpi_launcher.cmake")
include("${gemini_source}/test/hdf5_dll.cmake")
include("${gemini_source}/test/config.cmake")

set(WIN32 ${test_windows})
set(CMAKE_INSTALL_FULL_BINDIR "${CMAKE_CURRENT_BINARY_DIR}/local/bin")
set(MPIEXEC_EXECUTABLE "${CMAKE_COMMAND}")
set(MPIEXEC_NUMPROC_FLAG "-n")
set(MPIEXEC_MAX_NUMPROCS 2)
set(gemini3d_test_run_root "${CMAKE_CURRENT_BINARY_DIR}/native runs")
set(gemini3d_glow ON)
set(gemini3d_msis2 ON)

# Imported placeholders exercise the real CTest registration without native libraries.
foreach(target IN ITEMS gemini.bin gemini_c.bin gemini3d.run)
  add_executable(${target} IMPORTED)
  set_property(TARGET ${target} PROPERTY IMPORTED_LOCATION "${CMAKE_COMMAND}")
endforeach()

function(compare_gemini_output name outdir refdir)
endfunction()

include("${gemini_source}/test/sim.cmake")

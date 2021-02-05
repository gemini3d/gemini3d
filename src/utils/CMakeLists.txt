include(CheckSymbolExists)
include(CheckIncludeFile)

add_library(helpmod OBJECT help.f90)
target_link_libraries(helpmod PRIVATE config)
set_target_properties(helpmod PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(helpmod INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

#--------magnetic field calculation executable---------------
configure_file(magcalc_cli.in.f90 magcalc_cli.f90 @ONLY)
add_executable(magcalc.bin)
target_sources(magcalc.bin PRIVATE magcalc.f90 ${CMAKE_CURRENT_BINARY_DIR}/magcalc_cli.f90)
target_link_libraries(magcalc.bin PRIVATE helpmod helpmod config const grid io mesh mpimod pathlib reader sanity_check timeutils h5fortran::h5fortran)
set_target_properties(magcalc.bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

add_test(NAME unit:magcalc_exe_ok COMMAND $<TARGET_FILE:magcalc.bin>)
set_tests_properties(unit:magcalc_exe_ok PROPERTIES PASS_REGULAR_EXPRESSION "EOF: magcalc.bin")

# --- gemini3d.run Fortran front end
configure_file(gemini3d_run.in.f90 gemini3d_run.f90 @ONLY)
add_executable(gemini3d.run)
target_sources(gemini3d.run PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/gemini3d_run.f90 runner.f90 get_cpu.f90)
target_link_libraries(gemini3d.run PRIVATE const helpmod config pathlib timeutils autogrid reader h5fortran::h5fortran nc4fortran::nc4fortran)
set_target_properties(gemini3d.run PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set_target_properties(gemini3d.run PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_test(NAME unit:gemini3d_run_ok COMMAND $<TARGET_FILE:gemini3d.run>)
set_tests_properties(unit:gemini3d_run_ok PROPERTIES PASS_REGULAR_EXPRESSION "EOF: gemini3d.run")

function(cpu_count)
# Ref: https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/ProcessorCount.cmake
# https://gitlab.kitware.com/search?group_id=415&project_id=541&repository_ref=master&scope=blobs&search=GetNumberOfPhysicalCPU&snippets=false
# without hwloc, hyperthreading may give incorrect result
if(mpi)
  find_package(HWLOC)

  if(HWLOC_FOUND)
    # MSYS2 also has this
    target_sources(gemini3d.run PRIVATE get_cpu_hwloc.c)
    target_link_libraries(gemini3d.run PRIVATE HWLOC::HWLOC)
    return()
  endif()

  check_include_file(unistd.h HAVE_UNISTD_H)
  if(HAVE_UNISTD_H)
    check_symbol_exists(_SC_NPROCESSORS_ONLN unistd.h HAVE_GNU_COUNT)
  endif()

  if(HAVE_GNU_COUNT)
    # https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/sysconf.3.html
    target_sources(gemini3d.run PRIVATE get_cpu_sysconf.c)
  elseif(WIN32)
    target_sources(gemini3d.run PRIVATE get_cpu_win.c)
  else()
    target_sources(gemini3d.run PRIVATE dummy_cpu.c)
  endif()

else(mpi)
  target_sources(gemini3d.run PRIVATE dummy_cpu.c)
endif(mpi)

endfunction(cpu_count)

cpu_count()

# --- Gemini compare

add_library(assert assert.f90)
target_include_directories(assert INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(assert PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)


add_executable(gemini3d.compare)
target_sources(gemini3d.compare PRIVATE gemini_compare.f90 ${CMAKE_CURRENT_SOURCE_DIR}/../io/compare.f90)
target_link_libraries(gemini3d.compare PRIVATE assert config reader pathlib timeutils const nc4fortran::nc4fortran h5fortran::h5fortran)
target_include_directories(gemini3d.compare PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(gemini3d.compare PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

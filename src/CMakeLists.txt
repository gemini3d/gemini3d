function(set_targ_props)

foreach(t ${ARGV})
  target_include_directories(${t} INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(${t} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
endforeach()

endfunction(set_targ_props)


add_subdirectory(unit_tests)

# --- MSISE00 or MSIS 2.0
add_subdirectory(vendor/nrl_msis)

# --- HWM14
if(hwm14)
  add_subdirectory(vendor/nrl_hwm)
endif()

# --- just use one function from LAPACK95
add_subdirectory(vendor/lapack95)  # used by PDEparabolic

# --- GEMINI
add_subdirectory(ionization)
add_subdirectory(numerical)
# we use this even with MPI stub
add_subdirectory(mpimod)

add_subdirectory(io)
add_subdirectory(temporal)

add_library(collision OBJECT collisions/collisions.f90)
target_link_libraries(collision PRIVATE const config)

# --- other libs

add_subdirectory(neutral)
add_subdirectory(sources)

add_library(multifluid OBJECT multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE advec calculus collision const diffusion grid mesh ionization mpimod precipBCs sources timeutils config)

# --- Main GEMINI executable
configure_file(gemini_cli.in.f90 gemini_cli.f90 @ONLY)
add_executable(gemini.bin gemini.f90 ${CMAKE_CURRENT_BINARY_DIR}/gemini_cli.f90)
set_target_properties(gemini.bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_link_libraries(gemini.bin PRIVATE advec calculus config const collision diffusion fang grid io interp ionization gbsv mesh mpimod
  msis_ifc multifluid
  neutral pathlib PDEelliptic PDEparabolic potential precipBCs reader sanity_check sources temporal timeutils)
target_link_libraries(gemini.bin PRIVATE MUMPS::MUMPS MPI::MPI_Fortran MPI::MPI_C)
# Note: CMake 3.13 needs MPI linked here explicity to gemini.bin on any OS

add_test(NAME unit:gemini_exe_ok COMMAND $<TARGET_FILE:gemini.bin>)
set_tests_properties(unit:gemini_exe_ok PROPERTIES PASS_REGULAR_EXPRESSION "EOF: Gemini-3D")

#--------magnetic field calculation executable---------------
configure_file(magcalc_cli.in.f90 magcalc_cli.f90 @ONLY)
add_executable(magcalc.bin magcalc.f90 ${CMAKE_CURRENT_BINARY_DIR}/magcalc_cli.f90)
target_link_libraries(magcalc.bin PRIVATE config const grid io mesh mpimod pathlib reader sanity_check timeutils h5fortran::h5fortran)
set_target_properties(magcalc.bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

add_test(NAME unit:magcalc_exe_ok COMMAND $<TARGET_FILE:magcalc.bin>)
set_tests_properties(unit:magcalc_exe_ok PROPERTIES PASS_REGULAR_EXPRESSION "EOF: MAGCALC")

# --- gemini3d.run Fortran front end
add_executable(gemini3d.run gemini_runner.f90 mpimod/runner.f90 mpimod/get_cpu.f90)
target_link_libraries(gemini3d.run PRIVATE const config pathlib timeutils autogrid reader h5fortran::h5fortran nc4fortran::nc4fortran)
set_target_properties(gemini3d.run PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set_target_properties(gemini3d.run PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# Ref: https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/ProcessorCount.cmake
# https://gitlab.kitware.com/search?group_id=415&project_id=541&repository_ref=master&scope=blobs&search=GetNumberOfPhysicalCPU&snippets=false
# shortcomings: hyperthreading may give incorrect result
if(mpi)
  if(WIN32)
    # Windows
    target_sources(gemini3d.run PRIVATE mpimod/get_cpu_win.c)
  else()
    # MacOS and Linux
    find_path(hwloc_inc NAMES hwloc.h)
    find_library(hwloc NAMES hwloc)
    if(hwloc_inc AND hwloc)
      target_sources(gemini3d.run PRIVATE mpimod/get_cpu_hwloc.c)
      target_include_directories(gemini3d.run PRIVATE ${hwloc_inc})
      target_link_libraries(gemini3d.run PRIVATE ${hwloc})
    else()
      target_sources(gemini3d.run PRIVATE mpimod/get_cpu_sysconf.c)
    endif()
  endif()
else(mpi)
  target_sources(gemini3d.run PRIVATE mpimod/dummy_cpu.c)
endif(mpi)

# --- Gemini compare

if(${PROJECT_NAME}_BUILD_TESTING)
  add_executable(gemini3d.compare gemini_compare.f90 io/read_plasma_hdf5.f90)
  target_link_libraries(gemini3d.compare PRIVATE assert config reader pathlib timeutils const nc4fortran::nc4fortran h5fortran::h5fortran)
  target_include_directories(gemini3d.compare PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(gemini3d.compare PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
endif()

# --- apply properties

set_targ_props(autogrid calculus collision config const fang fang_run io ionization gbsv logging mpimod multifluid neutral pathlib precipBCs reader sanity_check sources temporal timeutils)

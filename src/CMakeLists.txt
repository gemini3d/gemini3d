add_subdirectory(unit_tests)

# --- MSISE00 or MSIS 2.0
add_subdirectory(vendor/nrl_msis)

# --- HWM14
add_subdirectory(vendor/nrl_hwm)

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
target_link_libraries(multifluid PRIVATE advec calculus collision const diffusion grid meshobj ionization mpimod precipBCs sources timeutils config)

add_subdirectory(utils)

# --- Main GEMINI3D library
configure_file(gemini_cli.in.f90 gemini_cli.f90 @ONLY)

add_library(gemini3d
libgemini.f90 gemini_init.f90 ${CMAKE_CURRENT_BINARY_DIR}/gemini_cli.f90
)
set_target_properties(gemini3d PROPERTIES
ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)
foreach(l
  advec calculus config const collision diffusion fang interp ionization gbsv
  meshobj meshobj_cart meshobj_dipole spherical newton geomagnetic
  msis_ifc multifluid
  neutral
  pathlib
  PDEparabolic precipBCs reader sanity_check errors sources temporal timeutils
  exe_frontend
  grid io PDEelliptic potential mpimod
  )
  target_link_libraries(gemini3d PUBLIC
  #$<TARGET_OBJECTS:${l}>
  ${l}
  )
endforeach()
target_link_libraries(gemini3d PRIVATE MPI::MPI_Fortran)

# --- Gemini.bin minimal main Fortran program
add_executable(gemini.bin gemini_main.f90
)
set_target_properties(gemini.bin PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)
target_link_libraries(gemini.bin PRIVATE
gemini3d
autogrid
hwloc_ifc hwloc_c
)

if(mpi)
  target_link_libraries(gemini.bin PRIVATE MUMPS::MUMPS SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)
else()
  target_link_libraries(gemini.bin PRIVATE MUMPS::MUMPS LAPACK::LAPACK MPI::MPI_Fortran)
endif()
if(MUMPS_HAVE_Scotch)
  target_link_libraries(gemini.bin PRIVATE Scotch::Scotch)
  target_link_libraries(gemini.bin PRIVATE METIS::METIS)
endif()
# same link logic applied to src/unit_tests/test_mumps

add_test(NAME unit:gemini_exe_ok COMMAND $<TARGET_FILE:gemini.bin>)
set_tests_properties(unit:gemini_exe_ok PROPERTIES
  PASS_REGULAR_EXPRESSION "EOF: gemini.bin"
  FIXTURES_SETUP gemini_exe_fxt
  FIXTURES_REQUIRED "mumps_fxt;hwloc_fxt"
  LABELS unit)

if(test_dll_path)
  set_tests_properties(unit:gemini_exe_ok PROPERTIES
    ENVIRONMENT "PATH=${test_dll_path}")
endif()

install(TARGETS gemini.bin
  EXPORT ${PROJECT_NAME}Targets)


# --- Gemini_c.bin
# This is an example of a C/C++ main program, instead of or in addition to the main Fortran gemini.bin
# this could be used as a starting point for C/C++ programs that use the GEMINI3D API

if(NOT mpi)
  message(STATUS "SKIP: GEMINI3D C program requires MPI")
  return()
endif()

add_executable(gemini_c.bin gemini_main.c)
set_target_properties(gemini_c.bin PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
LINKER_LANGUAGE C
)
# LINKER_LANGUAGE needed for Intel compilers
target_link_libraries(gemini_c.bin PRIVATE
gemini3d
autogrid
hwloc_ifc hwloc_c
MUMPS::MUMPS SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)

install(FILES gemini3d.h TYPE INCLUDE)

# NOTE: gemini_c.bin does not yet have unit testing, as we are deciding whether to make it a separate example
# or even in a separate project. See https://github.com/gemini3d/libgemini

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
add_subdirectory(inputdata)

# --- Main GEMINI3D library

configure_file(gemini_cli.in.f90 gemini_cli.f90 @ONLY)

add_library(gemini3d
libgemini.f90 gemini_init.f90 ${CMAKE_CURRENT_BINARY_DIR}/gemini_cli.f90
)
set_target_properties(gemini3d PROPERTIES
ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

set(gemini_internal_mods
advec
calculus config
const
collision
diffusion fang interp
ionization
inputdata precipdata efielddata
neutraldata neutraldata3D neutraldata2D neutraldata2Daxisymm neutraldata2Dcart
meshobj meshobj_cart meshobj_dipole spherical newton geomagnetic
multifluid
neutral
PDEparabolic gbsv
precipBCs reader sanity_check errors sources temporal timeutils
exe_frontend
grid io
PDEelliptic potential mumps_interface
mpimod
autogrid
msis_ifc msis00mod
hwm_ifc
hwloc_ifc hwloc_c
)
if(msis2)
  list(APPEND gemini_internal_mods msis2mod)
endif()

foreach(l ${gemini_internal_mods})

  install(TARGETS ${l}
  EXPORT ${PROJECT_NAME}-targets
  )

  target_sources(gemini3d PRIVATE $<TARGET_OBJECTS:${l}>)

  target_link_libraries(gemini3d
  # PRIVATE $<TARGET_OBJECTS:${l}>
  # $<BUILD_INTERFACE:${l}>
  $<INSTALL_INTERFACE:${l}>
  # ${l}
  )
endforeach()

# --- Gemini.bin minimal main Fortran program
add_executable(gemini.bin gemini_main.f90
)
set_target_properties(gemini.bin PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
RELWITHDEBINFO_POSTFIX .debug
DEBUG_POSTFIX .debug
)
target_link_libraries(gemini.bin PRIVATE
gemini3d
MUMPS::MUMPS
$<$<BOOL:${mpi}>:SCALAPACK::SCALAPACK>
LAPACK::LAPACK
nc4fortran::nc4fortran
h5fortran::h5fortran
$<$<BOOL:${HWLOC_FOUND}>:HWLOC::HWLOC>
MPI::MPI_Fortran
$<$<NOT:$<BOOL:${mpi}>>:MUMPS::MPISEQ>
"$<$<BOOL:${MUMPS_HAVE_Scotch}>:Scotch::Scotch;METIS::METIS>"
$<$<BOOL:${glow}>:glow::glow>
FortranFilesystem::filesystem
)

# same link logic applied to src/unit_tests/test_mumps

add_test(NAME unit:gemini_exe_ok COMMAND gemini.bin)
set_tests_properties(unit:gemini_exe_ok PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: gemini.bin"
FIXTURES_SETUP gemini_exe_fxt
FIXTURES_REQUIRED "mumps_fxt;hwloc_fxt"
LABELS unit
ENVIRONMENT $<$<BOOL:${test_dll_path}>:"PATH=${test_dll_path}">
)

install(TARGETS gemini3d gemini.bin
EXPORT ${PROJECT_NAME}-targets
)


# --- Gemini_c.bin
# This is an example of a C++ main program, instead of or in addition to the main Fortran gemini.bin
# this could be used as a starting point for C++ programs that use the GEMINI3D API

if(cpp_frontend)

include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/iniparser.cmake)
# for gemini_c C++ front end

add_executable(gemini_c.bin gemini_main.cpp pathlib.cpp)
set_target_properties(gemini_c.bin PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
#LINKER_LANGUAGE CXX
RELWITHDEBINFO_POSTFIX .debug
DEBUG_POSTFIX .debug
)
# LINKER_LANGUAGE needed for Intel compilers
target_link_libraries(gemini_c.bin PRIVATE
gemini3d
MUMPS::MUMPS
$<$<BOOL:${mpi}>:SCALAPACK::SCALAPACK>
LAPACK::LAPACK
nc4fortran::nc4fortran
h5fortran::h5fortran
$<$<BOOL:${HWLOC_FOUND}>:HWLOC::HWLOC>
MPI::MPI_Fortran
$<$<NOT:$<BOOL:${mpi}>>:MUMPS::MPISEQ>
"$<$<BOOL:${MUMPS_HAVE_Scotch}>:Scotch::Scotch;METIS::METIS>"
$<$<BOOL:${glow}>:glow::glow>
iniparser
)
target_compile_definitions(gemini_c.bin PRIVATE $<$<BOOL:${MSVC}>:_CRT_SECURE_NO_WARNINGS>)
target_compile_features(gemini_c.bin PRIVATE cxx_std_17)

# NOTE: gemini_c.bin does not yet have unit testing, as we are deciding whether to make it a separate example
# or even in a separate project. See https://github.com/gemini3d/libgemini

install(TARGETS gemini_c.bin
EXPORT ${PROJECT_NAME}-targets
)
install(FILES gemini3d.h TYPE INCLUDE)

endif(cpp_frontend)

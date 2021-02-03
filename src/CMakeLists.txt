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

# --- apply properties

set_targ_props(autogrid calculus collision config const fang fang_run io ionization gbsv logging mpimod multifluid neutral pathlib precipBCs reader sanity_check sources temporal timeutils)

# --- extra programs

add_subdirectory(utils)

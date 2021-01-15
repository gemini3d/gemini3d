function(set_targ_props)

foreach(t ${ARGV})
  target_include_directories(${t} INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(${t} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
endforeach()

endfunction(set_targ_props)

# --- MSISE00 or MSIS 2.0
add_subdirectory(vendor/msis00)
if(msis20)
  # provides gtd8() subroutine
  include(${PROJECT_SOURCE_DIR}/cmake/msis2.cmake)
endif()

# --- HWM14
if(hwm14)
  include(${PROJECT_SOURCE_DIR}/cmake/hwm14.cmake)
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
if(mpi AND ${PROJECT_NAME}_BUILD_TESTING)
  # keep this if(mpi) because these tests truly use MPI -- don't work with MPI stub
  add_subdirectory(unit_tests)
endif()

add_subdirectory(neutral)
add_subdirectory(sources)

add_library(multifluid OBJECT multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE advec calculus collision const diffusion grid mesh ionization mpimod precipBCs sources timeutils config)

# --- Main GEMINI executable
configure_file(gemini_cli.in.f90 gemini_cli.f90 @ONLY)
add_executable(gemini.bin gemini.f90 ${CMAKE_CURRENT_BINARY_DIR}/gemini_cli.f90)
set_target_properties(gemini.bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_link_libraries(gemini.bin PRIVATE advec calculus config const collision diffusion fang grid io interp ionization gbsv mesh mpimod multifluid
  neutral pathlib PDEelliptic PDEparabolic potential precipBCs reader sanity_check sources temporal timeutils)
target_link_libraries(gemini.bin PRIVATE MUMPS::MUMPS)

add_test(NAME unit:gemini_exe_ok COMMAND $<TARGET_FILE:gemini.bin>)
set_tests_properties(unit:gemini_exe_ok PROPERTIES PASS_REGULAR_EXPRESSION "EOF: Gemini-3D")

#--------magnetic field calculation executable---------------
configure_file(magcalc_cli.in.f90 magcalc_cli.f90 @ONLY)
add_executable(magcalc.bin magcalc.f90 ${CMAKE_CURRENT_BINARY_DIR}/magcalc_cli.f90)
target_link_libraries(magcalc.bin PRIVATE config const grid io mesh mpimod pathlib reader sanity_check timeutils h5fortran::h5fortran)
set_target_properties(magcalc.bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

add_test(NAME unit:magcalc_exe_ok COMMAND $<TARGET_FILE:magcalc.bin>)
set_tests_properties(unit:magcalc_exe_ok PROPERTIES PASS_REGULAR_EXPRESSION "EOF: MAGCALC")

set_targ_props(calculus collision config const fang fang_run io ionization gbsv logging mpimod multifluid neutral pathlib precipBCs reader sanity_check sources temporal timeutils)

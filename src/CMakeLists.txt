# --- simple tests to help catch broken libs
add_subdirectory(tests)

# --- vendor libraries
add_subdirectory(vendor/msis00)
add_subdirectory(vendor/lapack95)  # used by PDEparabolic

# --- GEMINI
add_subdirectory(ionization)
add_subdirectory(numerical)
add_subdirectory(io)
add_subdirectory(temporal)
add_subdirectory(neutral)

add_library(collision OBJECT collisions/collisions.f90)
target_link_libraries(collision PRIVATE const)

add_subdirectory(sources)

# --- other libs
add_library(multifluid OBJECT multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE advec calculus collision config const diffusion grid mesh ionization mpimod precipBCs sources timeutils)

# --- Main GEMINI executable
add_executable(gemini.bin gemini.f90)
target_link_libraries(gemini.bin PRIVATE advec calculus config const collision diffusion fang grid io interp ionization gbsv mesh mpimod multifluid
  neutral pathlib PDEelliptic PDEparabolic potential precipBCs reader sources temporal timeutils
  ${LAPACK_LIBRARIES})  # necessary for self-compiled
set_target_properties(gemini.bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#--------magnetic field calculation executable---------------
add_executable(magcalc.bin magcalc.f90)
target_link_libraries(magcalc.bin PRIVATE config const grid io mesh mpimod pathlib reader timeutils)
set_target_properties(magcalc.bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# --- set common target properties

foreach(t collision config fang fang_run io ionization gbsv logging multifluid neutral pathlib precipBCs reader sources temporal timeutils)
  target_include_directories(${t} INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(${t} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
endforeach()

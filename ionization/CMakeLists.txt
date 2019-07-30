cmake_policy(SET CMP0079 NEW)
cmake_policy(SET CMP0076 NEW)

target_link_libraries(ionization PRIVATE neutral const mpimod)
set_target_properties(ionization PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
if(USEGLOW)
  include(${PROJECT_SOURCE_DIR}/cmake/glow.cmake)

  target_sources(ionization PRIVATE glow_run.F90)
  target_link_libraries(ionization PRIVATE cglow)
  target_include_directories(ionization PRIVATE ${ncarglow_BINARY_DIR}/src)
  target_compile_definitions(ionization PRIVATE
    DATADIR="build/_deps/ncarglow-src/data/")
else()
  target_sources(ionization PRIVATE glow_dummy.f90)
endif(USEGLOW)


add_library(precipBCs boundary_conditions/precipBCs_mod.f90)
target_link_libraries(precipBCs PRIVATE const grid io interp temporal timeutils)
target_compile_options(precipBCs PRIVATE ${FFLAGS})
set_target_properties(precipBCs PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
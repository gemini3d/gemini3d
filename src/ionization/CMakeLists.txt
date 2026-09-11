set_property(DIRECTORY PROPERTY LABELS ionization)

add_library(fang OBJECT fang.f90)
target_link_libraries(fang PRIVATE const)

add_library(fang_run OBJECT fang_run.f90)
target_link_libraries(fang_run PRIVATE fang msis::msis const)

add_library(ionization OBJECT ionization.f90
$<TARGET_OBJECTS:grid>
)

if(gemini3d_glow)
  # must have trailing slash in path
  if(NOT DEFINED glow_data_dir)
    message(FATAL_ERROR "[GEMINI3D] GLOW data directory variable is not defined: glow_data_dir")
  endif()

  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(glow_data_dir_free_flag "-ffree-line-length-none")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
    # Flang does not have a free-form column limit https://github.com/llvm/llvm-project/issues/167912#issuecomment-3529004167
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    # this diagnostic occurs if compiler option "-stand" was used
    set(glow_data_dir_free_flag "-diag-disable:5268")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(glow_data_dir_free_flag "-Mextend")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    # already allows no limit for free form
  endif()

  configure_file(glow_run.in.f90 glow_run.f90 @ONLY)
  set_property(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/glow_run.f90 APPEND PROPERTY COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:Fortran>:${glow_data_dir_free_flag}>)

  target_sources(ionization PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/glow_run.f90)
else()
  target_sources(ionization PRIVATE glow_dummy.f90)
endif()

target_link_libraries(ionization PRIVATE const meshobj timeutils fang
$<$<BOOL:${gemini3d_glow}>:glow::glow>
)

add_library(ionization_mpi OBJECT ionization_mpi.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(ionization_mpi PRIVATE ionization const neutral meshobj timeutils fang
mpimod
MPI::MPI_Fortran
)

add_library(precipBCs OBJECT boundary_conditions/precipBCs_mod.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(precipBCs PRIVATE reader const meshobj interp timeutils gemini3d_config
inputdata precipdata
MPI::MPI_Fortran
)

add_library(solfluxBCs OBJECT boundary_conditions/solfluxBCs_mod.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(solfluxBCs PRIVATE reader const meshobj interp timeutils gemini3d_config
inputdata solfluxdata
MPI::MPI_Fortran
)

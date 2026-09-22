add_executable(audit_cartmesh_metrics cartmesh_metrics.f90)
target_link_libraries(audit_cartmesh_metrics PRIVATE meshobj_cart const gemini3d)
intel_fortran_main_linker(audit_cartmesh_metrics)
add_test(NAME audit:cartmesh_metrics COMMAND audit_cartmesh_metrics)
set_tests_properties(audit:cartmesh_metrics PROPERTIES LABELS "unit;audit;numerics")

add_executable(audit_halo_end halo_end.f90)
target_link_libraries(audit_halo_end PRIVATE mpimod const gemini3d MPI::MPI_Fortran)
intel_fortran_main_linker(audit_halo_end)

# magcalc keeps these procedures internal to its main program. Extract them
# without edits so the regression exercises production code, not a test copy.
set(_magcalc_source ${PROJECT_SOURCE_DIR}/src/utils/magcalc.f90)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_magcalc_source})
file(READ ${_magcalc_source} _magcalc)
string(REGEX MATCH "  subroutine fixJ\\(J1,J2,J3\\).*end function integrate2D" MAGCALC_KERNELS "${_magcalc}")
string(REGEX MATCH "      Rcubed\\(:,:,:\\)=Rx\\*\\*2\\+Ry\\*\\*2.*Bphi\\(ipoints\\) = integrate2D\\(integrand,integrandend\\)"
  MAGCALC_2D "${_magcalc}")
if(NOT MAGCALC_KERNELS OR NOT MAGCALC_2D)
  message(FATAL_ERROR "Could not extract production magcalc kernels for audit regression")
endif()
configure_file(magcalc_kernels.f90.in magcalc_kernels.f90 @ONLY)
add_executable(audit_magcalc_kernels ${CMAKE_CURRENT_BINARY_DIR}/magcalc_kernels.f90)
target_link_libraries(audit_magcalc_kernels PRIVATE mpimod meshobj_cart const gemini3d MPI::MPI_Fortran)
intel_fortran_main_linker(audit_magcalc_kernels)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(audit_magcalc_kernels PRIVATE -ffpe-trap=invalid,zero,overflow)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel" AND NOT WIN32)
  target_compile_options(audit_magcalc_kernels PRIVATE -fpe0)
endif()
foreach(ranks IN ITEMS 1 2 4)
  foreach(probe IN ITEMS halo_end magcalc_kernels)
    add_test(NAME audit:${probe}_${ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:audit_${probe}> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(audit:${probe}_${ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${ranks} LABELS "unit;audit;numerics;mpi")
  endforeach()
endforeach()

add_executable(audit_diffusion_singular diffusion_singular.f90)
target_link_libraries(audit_diffusion_singular PRIVATE PDEparabolic const MPI::MPI_Fortran)
if(Python_Interpreter_FOUND)
  foreach(stage IN ITEMS euler tr bdf2)
    add_test(NAME audit:diffusion_singular_${stage}
      COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_diffusion_failure.py
        $<TARGET_FILE:audit_diffusion_singular> ${stage})
    set_tests_properties(audit:diffusion_singular_${stage} PROPERTIES TIMEOUT 30 LABELS "unit;audit;numerics")
  endforeach()
endif()

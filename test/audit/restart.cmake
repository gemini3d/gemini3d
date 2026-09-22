# Focused synthetic-state checks; these do not qualify research trajectories.
add_executable(qualification_restart_runtime restart_probe.f90)
target_link_libraries(qualification_restart_runtime PRIVATE restart_runtime atomic_file const
  gemini3d_config timeutils mpimod gemini3d h5fortran::h5fortran MPI::MPI_Fortran)
intel_fortran_main_linker(qualification_restart_runtime)
if(NUMPY_FOUND AND H5PY_FOUND)
  set(restart_mpi_flags)
  foreach(flag IN LISTS MPIEXEC_PREFLAGS)
    list(APPEND restart_mpi_flags "--mpi-pref=${flag}")
  endforeach()
  foreach(flag IN LISTS MPIEXEC_POSTFLAGS)
    list(APPEND restart_mpi_flags "--mpi-post=${flag}")
  endforeach()
  foreach(layout IN ITEMS 1x2 2x2)
    if(layout STREQUAL "2x2")
      set(ranks 4)
    else()
      set(ranks 2)
    endif()
    set(layout_restart_mpi_flags ${restart_mpi_flags})
    if(ranks GREATER 2 AND MPI_C_LIBRARY_VERSION_STRING MATCHES "Open[ ]?MPI")
      list(APPEND layout_restart_mpi_flags "--mpi-pref=--map-by" "--mpi-pref=:OVERSUBSCRIBE")
    endif()
    add_test(NAME qualification:restart_runtime_${layout} COMMAND ${Python_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/test/qualification/test_restart_runtime.py
      --exe $<TARGET_FILE:qualification_restart_runtime>
      --mpiexec ${MPIEXEC_EXECUTABLE} --layout ${layout}
      --numproc-flag=${MPIEXEC_NUMPROC_FLAG} ${layout_restart_mpi_flags}
      --work ${CMAKE_CURRENT_BINARY_DIR}/restart-runtime-${layout})
    set_tests_properties(qualification:restart_runtime_${layout} PROPERTIES
      LABELS "unit;qualification;mpi" TIMEOUT 180 PROCESSORS ${ranks}
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1;PRTE_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

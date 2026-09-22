add_executable(audit_neutral_background neutral_background_test.f90)
target_link_libraries(audit_neutral_background PRIVATE neutral_background neutral neutraldataBG meshobj_cart grid gemini3d)
intel_fortran_main_linker(audit_neutral_background)
foreach(mode IN ITEMS covered endpoint singleton descending_covered ascending descending closed closed_covered zero belowground)
  add_test(NAME audit:neutral_background_${mode} COMMAND audit_neutral_background ${mode})
  set_tests_properties(audit:neutral_background_${mode} PROPERTIES LABELS "unit;audit")
endforeach()
if(Python_Interpreter_FOUND)
  foreach(mode IN ITEMS lower_coverage upper_coverage zero_denominator negative_density underground temperature singleton_uncovered)
    add_test(NAME audit:neutral_background_reject_${mode} COMMAND ${Python_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/neutral_rejection.py $<TARGET_FILE:audit_neutral_background> ${mode})
    set_tests_properties(audit:neutral_background_reject_${mode} PROPERTIES LABELS "unit;audit")
  endforeach()
endif()

add_executable(audit_input_lifecycle input_lifecycle_test.f90)
target_link_libraries(audit_input_lifecycle PRIVATE inputdata precipdata solfluxdata efielddata neutraldataBG
  neutraldata neutraldata2D neutraldata2Dcart neutraldata2Daxisymm neutraldata3D neutraldata3D_mpi
  neutraldata3D_geom_mpi neutraldata3D_geog_mpi neutraldata3D_fclaw neutraldata3D_fclaw_axisymm
  neutraldata3D_fclaw_3Dx meshobj_cart meshobj_dipole gemini3d)
intel_fortran_main_linker(audit_input_lifecycle)
add_test(NAME audit:input_lifecycle COMMAND audit_input_lifecycle)
set_tests_properties(audit:input_lifecycle PROPERTIES LABELS "unit;audit")

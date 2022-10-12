set_property(DIRECTORY PROPERTY  LABELS ionization)

add_executable(test_fang test_fang.f90)
target_link_libraries(test_fang PRIVATE fang fang_run assert const msis::msis_ifc)

add_test(NAME fangIonize COMMAND test_fang)
set_property(TEST fangIonize PROPERTY LABELS unit)
dll_test_path(msis::msis_ifc fangIonize)

if(H5PY_FOUND)
  add_test(NAME fangIonize:python
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_fang.py $<TARGET_FILE:test_fang>
  )
  set_property(TEST fangIonize:python PROPERTY LABELS "unit;python")
endif()

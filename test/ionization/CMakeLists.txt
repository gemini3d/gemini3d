set_directory_properties(PROPERTIES LABELS ionization)

add_executable(test_fang test_fang.f90)
target_link_libraries(test_fang PRIVATE fang fang_run assert const msis::msis_ifc)

add_test(NAME fangIonize COMMAND test_fang)
set_tests_properties(fangIonize PROPERTIES
TIMEOUT 10
LABELS unit
)
dll_test_path(msis::msis_ifc fangIonize)

if(python)
  add_test(NAME fangIonize:python
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_fang.py $<TARGET_FILE:test_fang>
  )
  set_tests_properties(fangIonize:python PROPERTIES
  TIMEOUT 30
  LABELS "unit;python"
  )
endif()

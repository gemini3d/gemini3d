set_property(DIRECTORY PROPERTY  LABELS ionization)

add_executable(test_fang test_fang.f90)
target_link_libraries(test_fang PRIVATE fang fang_run assert const msis::msis_ifc)

add_test(NAME fangIonize:MSISE00 COMMAND test_fang 0)
add_test(NAME fangIonize:MSIS21 COMMAND test_fang 21)

set_tests_properties(fangIonize:MSISE00 fangIonize:MSIS21 PROPERTIES
TIMEOUT 10
LABELS unit
WORKING_DIRECTORY $<TARGET_FILE_DIR:gemini.bin>
)
dll_test_path(msis::msis_ifc "fangIonize:MSISE00;fangIonize:MSIS21")

add_test(NAME fangIonize:python
COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_fang.py $<TARGET_FILE:test_fang>
)
set_tests_properties(fangIonize:python PROPERTIES
TIMEOUT 30
LABELS "unit;python"
DISABLED $<NOT:$<BOOL:${H5PY_FOUND}>>
)

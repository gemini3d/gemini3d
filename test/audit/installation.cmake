if(Python_Interpreter_FOUND AND Python_VERSION VERSION_GREATER_EQUAL 3.11)
  add_test(NAME qualification:local_environment
    COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/qualification/test_local_environment.py)
  set_tests_properties(qualification:local_environment PROPERTIES LABELS "unit;qualification;installation")
endif()

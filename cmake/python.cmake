if(NOT DEFINED PythonOK)

find_package(Python3 COMPONENTS Interpreter)
# Python3::Interpreter did NOT work
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import gemini; print(gemini.__file__)"
  ERROR_QUIET OUTPUT_QUIET
  RESULT_VARIABLE ret
  TIMEOUT 15)
if(ret EQUAL 0)
  set(PythonOK true CACHE BOOL "PyGemini is present.")
else()
  message(WARNING "Need to setup PyGemini by 'pip install -e gemini'   error: ${ret}")
  set(PythonOK false CACHE BOOL "PyGemini is NOT present.")
endif()

endif()
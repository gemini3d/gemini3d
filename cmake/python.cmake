
find_package(Python3 COMPONENTS Interpreter)

# Python3::Interpreter does NOT work, use ${Python3_EXECUTABLE}
if(NOT DEFINED python_ok)

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import gemini3d; print(gemini3d.__file__)"
  OUTPUT_VARIABLE _pygemloc
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE python_ok
  TIMEOUT 15)

if(python_ok EQUAL 0)
  message(STATUS "PyGemini: ${_pygemloc}")
  set(python_ok true CACHE BOOL "PyGemini OK")
else()
  set(python_ok false)
  message(STATUS "MISSING: PyGemini -> install by 'python3 setup.py --user develop'")
endif()

endif(NOT DEFINED python_ok)
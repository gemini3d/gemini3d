
find_package(Python3 COMPONENTS Interpreter)

# Python3::Interpreter does NOT work, use ${Python3_EXECUTABLE}
if(NOT DEFINED python_disabled)

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import gemini3d; print(gemini3d.__file__)"
  OUTPUT_VARIABLE _pygemloc
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE python_disabled
  TIMEOUT 15)

if(python_disabled EQUAL 0)
  message(STATUS "PyGemini: ${_pygemloc}")
  set(python_disabled $${python_disabled} CACHE STRING "PyGemini OK")
else()
  message(STATUS "MISSING: PyGemini -> install by 'python3 setup.py --user develop'")
endif()

endif(NOT DEFINED python_disabled)
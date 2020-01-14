
find_package(Python3 COMPONENTS Interpreter)
# keep find_package outside if() statement or Python3_EXECUTABLE will intermittently be missing

# Python3::Interpreter did NOT work
if(NOT DEFINED PythonOK)

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import gemini; print(gemini.__file__)"
  RESULT_VARIABLE ret
  TIMEOUT 15)
if(ret EQUAL 0)
  set(PythonOK true CACHE BOOL "PyGemini is present.")
else()
  message(STATUS "MISSING: PyGemini -> install by 'python3 setup.py --user develop'")
  set(PythonOK false CACHE BOOL "PyGemini is NOT present.")
endif(ret EQUAL 0)

endif(NOT DEFINED PythonOK)
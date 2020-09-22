if(NOT DEFINED Matlab_FOUND)

find_package(Matlab COMPONENTS MAIN_PROGRAM)

execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch "setup; gemini3d.tests.test_unit"
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  TIMEOUT 120
  RESULT_VARIABLE _ok
  )

if(_ok EQUAL 0)
  message(STATUS "MatGemini working")
  set(Matlab_FOUND true CACHE BOOL "MatGemini not working")
else()
  message(STATUS "Could not use MatGemini")
  set(Matlab_FOUND false CACHE BOOL "MatGemini not working" FORCE)
endif()

endif(NOT DEFINED Matlab_FOUND)

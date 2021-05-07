# --- find generator
function(find_generator)

find_program(ninja NAMES ninja ninja-build samu)

if(ninja)
  execute_process(COMMAND ${ninja} --version
    OUTPUT_VARIABLE ninja_version OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE err
    TIMEOUT 5)
  if(err EQUAL 0 AND ninja_version VERSION_GREATER_EQUAL 1.10)
    set(CTEST_CMAKE_GENERATOR Ninja)
  endif()
endif(ninja)

if(NOT DEFINED CTEST_CMAKE_GENERATOR)
  set(CTEST_BUILD_FLAGS -j)  # not --parallel as this goes to generator directly
  if(WIN32)
    set(CTEST_CMAKE_GENERATOR "MinGW Makefiles")
  else()
    set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  endif()
endif()

set(CTEST_CMAKE_GENERATOR ${CTEST_CMAKE_GENERATOR} PARENT_SCOPE)

endfunction(find_generator)

if(NOT DEFINED CTEST_CMAKE_GENERATOR)
  if(DEFINED ENV{CMAKE_GENERATOR})
    set(CTEST_CMAKE_GENERATOR $ENV{CMAKE_GENERATOR})
  else()
    find_generator()
  endif()
endif()

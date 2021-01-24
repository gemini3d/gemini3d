# this file defines simulation tests.
# The names of these must match those in cmake/gemini3d_url.ini

set(sim_timeout 1800)  # seconds to allow simulation tests to run

# --- get metadata

if(download)
  if(CMAKE_VERSION VERSION_LESS 3.19)
    message(FATAL_ERROR "Reference data download requires CMake >= 3.19")
  endif()

  file(READ ${CMAKE_CURRENT_LIST_DIR}/gemini3d_url.json _refj)
endif()

# --- setup tests

set(_tests 2dns_fang 2dew_fang 3d_fang)

if(glow)
  list(APPEND _tests 2dns_glow 2dew_glow 3d_glow)
endif(glow)

foreach(_s ${_tests})
  if(download)

    string(JSON _url GET ${_refj} ${_s} url)

    message(STATUS "sync ${_s}")
    FetchContent_Declare(download_${_s}
      URL ${_url}
      SOURCE_DIR ${PROJECT_SOURCE_DIR}/tests/data/test${_s}
      SOURCE_SUBDIR test${_s}
      TLS_VERIFY ON
      UPDATE_DISCONNECTED ON)
    FetchContent_MakeAvailable(download_${_s})

  endif(download)

  setup_gemini_test(${_s} ${sim_timeout})

endforeach()

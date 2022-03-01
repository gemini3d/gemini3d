function(dll_test_path libs test_names)
# if shared lib on Windows, need DLL on PATH

if(NOT WIN32)
  return()
endif()

if(CMAKE_VERSION VERSION_LESS 3.22)
  message(VERBOSE "CMake ${CMAKE_VERSION} < 3.22: cannot apply ENVIRONMENT_MODIFICATION to ${test_names}")
  return()
endif()

set(dll_mod)

foreach(lib IN LISTS libs)

  get_target_property(ttype ${lib} TYPE)
  if(NOT ttype STREQUAL SHARED_LIBRARY)
    message(DEBUG "${lib} is not a shared library, no need for ENVIRONMENT_MODIFICATION for ${test_names}")
    continue()
  endif()

  get_target_property(conf ${lib} IMPORTED_CONFIGURATIONS)
  if(conf)
    list(GET conf 0 conf)
    # assume first configuration is desired

    get_target_property(loc ${lib} IMPORTED_LOCATION_${conf})
    if(NOT loc)
      message(VERBOSE "did not find imported location for ${lib}")
      continue()
    endif()

    cmake_path(GET loc PARENT_PATH loc)

  else()

    set(loc $<TARGET_FILE_DIR:${lib}>)

  endif()

  list(APPEND dll_mod "PATH=path_list_append:${loc}")

endforeach()

set_tests_properties(${test_names} PROPERTIES
ENVIRONMENT_MODIFICATION "${dll_mod}"
)


endfunction(dll_test_path)

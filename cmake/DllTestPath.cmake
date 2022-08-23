function(dll_test_path libs test_names)
# if shared lib on Windows, need DLL on PATH

if(NOT WIN32)
  return()
endif()

if(CMAKE_VERSION VERSION_LESS 3.22)
  message(VERBOSE " CMake ${CMAKE_VERSION} < 3.22: cannot apply ENVIRONMENT_MODIFICATION to ${test_names}")
  return()
endif()

set(dll_mod)

foreach(lib IN LISTS libs)

  get_target_property(ttype ${lib} TYPE)
  if(ttype STREQUAL "STATIC_LIBRARY")
    message(DEBUG "${lib} is ${ttype}. No need for ENVIRONMENT_MODIFICATION for ${test_names}")
    continue()
  endif()

  get_target_property(imploc ${lib} IMPORTED_LOCATION_RELEASE)
  get_target_property(intloc ${lib} INTERFACE_LINK_LIBRARIES)
  if(imploc)
    foreach(l IN LISTS imploc)
      cmake_path(GET l PARENT_PATH loc)
      if(IS_DIRECTORY ${loc})
        list(APPEND dll_mod "PATH=path_list_append:${loc}")
        cmake_path(SET d NORMALIZE ${loc}/../bin)
        # can't check bin/stem.dll as some libs add arbitrary stuff to stem
        if(IS_DIRECTORY ${d})
          list(APPEND dll_mod "PATH=path_list_append:${d}")
        endif()
      endif()
    endforeach()
  elseif(intloc)
    foreach(l IN LISTS intloc)
      cmake_path(GET l PARENT_PATH loc)
      if(IS_DIRECTORY ${loc})
        list(APPEND dll_mod "PATH=path_list_append:${loc}")
        cmake_path(SET d NORMALIZE ${loc}/../bin)
        if(IS_DIRECTORY ${d})
          list(APPEND dll_mod "PATH=path_list_append:${d}")
        endif()
      endif()
    endforeach()
  elseif(EXISTS ${lib})
    list(APPEND dll_mod "PATH=path_list_append:$<TARGET_FILE_DIR:${lib}>")
  else()
    message(DEBUG "did not find library for ${lib} for ${test_names}")
  endif()

endforeach()

list(REMOVE_DUPLICATES dll_mod)

if(dll_mod)
  message(VERBOSE "environment_modification ${dll_mod} for ${test_names}")

  set_tests_properties(${test_names} PROPERTIES
  ENVIRONMENT_MODIFICATION "${dll_mod}"
  )
else()
  message(VERBOSE "no environment_modification for ${test_names}")
endif()


endfunction(dll_test_path)

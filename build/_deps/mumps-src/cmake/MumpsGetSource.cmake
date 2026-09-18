# retrieve MUMPS source files from JSON

include_guard()

function(mumps_get_src out key1 key2)
  string(JSON L LENGTH "${source_json}" "mumps_sources" "${key1}" "${key2}")

  set(_files)
  if(L GREATER 0)
    math(EXPR L "${L} - 1")
    foreach(_idx RANGE ${L})
      string(JSON _src_file GET "${source_json}" "mumps_sources" "${key1}" "${key2}" ${_idx})
      list(APPEND _files ${_src_file})
    endforeach()
  endif()

  set(${out} "${_files}" PARENT_SCOPE)
endfunction()

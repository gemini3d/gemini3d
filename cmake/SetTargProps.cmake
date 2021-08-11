function(set_targ_props)

foreach(t ${ARGV})

  # get_target_property(_bin ${t} BINARY_DIR)
  set(_bin ${PROJECT_BINARY_DIR})

  target_include_directories(${t} INTERFACE
  $<BUILD_INTERFACE:${_bin}/include>
  $<INSTALL_INTERFACE:include>)
  set_target_properties(${t} PROPERTIES Fortran_MODULE_DIRECTORY ${_bin}/include)
endforeach()

endfunction(set_targ_props)

function(c_check)

set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

if(ffilesystem_trace)
check_symbol_exists(__has_include "" c23_has_include)
endif()

check_symbol_exists(__has_c_attribute "" c23_has_c_attribute)
if(c23_has_c_attribute)
  check_source_compiles(C
  "#if !__has_c_attribute(maybe_unused)
  #error \"no maybe_unused\"
  #endif"
  c23_maybe_unused)

  # [[reproducible]] [[unsequenced]] support by compilers:
  # https://en.cppreference.com/w/c/compiler_support/23

  # check_source_compiles(C
  # "#if !__has_c_attribute(reproducible)
  # #error \"no reproducible\"
  # #endif"
  # c23_reproducible)

  # check_source_compiles(C
  # "#if !__has_c_attribute(unsequenced)
  # #error \"no unsequenced\"
  # #endif"
  # c23_unsequenced)
endif()

endfunction()

# Opt-in developer instrumentation; never silently enable in a production build.
option(gemini3d_sanitizers "GNU AddressSanitizer and UndefinedBehaviorSanitizer instrumentation" OFF)
if(gemini3d_sanitizers)
  if(NOT CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
     NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
     NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "This sanitizer profile is qualified only with the GNU C/C++/Fortran toolchain")
  endif()
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

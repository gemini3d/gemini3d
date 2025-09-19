
# flags we don't want leaking into Git submodules to avoid excessive warnings on projects we don't control.
set(${PROJECT_NAME}_flags
$<$<AND:$<CONFIG:Debug,RelWithDebInfo>,$<COMPILE_LANGUAGE:Fortran>>:-warn>
$<$<AND:$<CONFIG:Debug,RelWithDebInfo>,$<COMPILE_LANGUAGE:C,CXX>>:-Wall>
$<$<COMPILE_LANGUAGE:Fortran>:-traceback>
)

# this flag needs to be applied EVERYWHERE incl. submodule projects
# or runtime error / weird behavior with non-standard C_BOOL values.
# -standard-semantics is no good because it breaks linkage within oneAPI itself e.g. oneMPI library!
if(NOT WIN32)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-fpscomp;logicals>")
endif()

# add_compile_options($<$<CONFIG:Debug>:-Rno-debug-disables-optimization>)
# remark #5415: Feature not yet implemented: Some 'check' options temporarily disabled.
# warning #10182: disabling optimization; runtime debug checks enabled
# remark #7712: This variable has not been used.  (-Wunused-dummy-argument)


if(NOT WIN32)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    list(APPEND ${PROJECT_NAME}_flags -fiopenmp)
  endif()
  # -fiopenmp:
  # undefined reference to `omp_get_max_threads'
  # undefined reference to `__kmpc_global_thread_num' and more similar
  list(APPEND ${PROJECT_NAME}_flags -qopenmp)
  # undefined reference to `__kmpc_begin'
endif()

list(APPEND ${PROJECT_NAME}_flags $<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays>)
# heap-arrays: avoid stack overflow, for both unit tests and actual simulations
# it's needed on Linux and Windows
# https://www.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/compiler-reference/compiler-options/advanced-optimization-options/heap-arrays.html


# --- IMPORTANT: bounds checking
# add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-check>")
# -check is an alias for -check all. However, MUMPS trips on -check, so we have to use a less stringent check.
list(APPEND ${PROJECT_NAME}_flags "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-CB>")
# -CB is an alias for -check bounds.
# --- IMPORTANT

list(APPEND ${PROJECT_NAME}_flags "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-debug>")
# -debug is an alias for -debug all
# -fpe0 causes MUMPS failures (internal to MUMPS)

# Fortran 2018 standard too many false warnings

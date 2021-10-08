# NOTE: don't use -march=native as GCC doesn't support all CPU arches with that option.
add_compile_options(-mtune=native)

add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fimplicit-none>)
if(dev)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Wall;-Wextra>")
  # -Wpedantic makes too many false positives
else(dev)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument;-Wno-unused-variable;-Wno-unused-function>")
endif(dev)

# avoid backtrace that's unusable without -g
add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-fno-backtrace>)

# Wdo-subscript is known to warn on obvious non-problems
check_compiler_flag(Fortran -Wdo-subscript dosubflag)
if(dosubflag)
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Wno-do-subscript>)
endif()

add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-Werror=array-bounds;-fcheck=all>")
  # add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-ffpe-trap=invalid,zero,overflow>")#,underflow)

# Fortran 2018 standard flag is buggy at least through GCC 9, causing fake "Common block" warnings

if(CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 9.3.0)
  # makes a lot of spurious warnngs on alloctable scalar character
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Wno-maybe-uninitialized>)
endif()

check_compiler_flag(Fortran -fallow-argument-mismatch allow_mismatch_args)
if(allow_mismatch_args)
  set(gfortran_opts -fallow-argument-mismatch)
endif()

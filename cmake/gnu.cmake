include(CheckCompilerFlag)
# NOTE: don't use -march=native as GCC doesn't support all CPU arches with that option.
# add_compile_options(-mtune=native)

# flags we don't want leaking into Git submodules to avoid excessive warnings on projects we don't control.
set(${PROJECT_NAME}_flags
$<$<CONFIG:Debug,RelWithDebInfo>:-Wall>
$<$<COMPILE_LANGUAGE:Fortran>:-fimplicit-none>
$<$<COMPILE_LANGUAGE:Fortran>:-Wno-uninitialized>
$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<VERSION_LESS:$<Fortran_COMPILER_VERSION>,10>>:-Wno-conversion>
$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_VERSION:9.3.0>>:-Wno-maybe-uninitialized>
)
# lot of spurious warnings on allocatable scalar character -Wno-maybe-uninitialized

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
   CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0" AND
   CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.2"
   )
  list(APPEND ${PROJECT_NAME}_flags "$<$<COMPILE_LANGUAGE:Fortran>:-Wno-external-argument-mismatch>")
endif()
# workaround for src/numerical/coord/newton.f90 build fail despite the code being correct.
# https://www.scivision.dev/gfortran-15-external-argument-mismatch/

# --- IMPORTANT
list(APPEND ${PROJECT_NAME}_flags "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-Werror=array-bounds$<SEMICOLON>-fcheck=all>")
# --- IMPORTANT: options help trap array indexing/bounds errors at runtime

# avoid backtrace that's unusable without -g
list(APPEND ${PROJECT_NAME}_flags "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-fno-backtrace>")

# Wdo-subscript is known to warn on obvious non-problems
check_compiler_flag(Fortran -Wdo-subscript dosubflag)
if(dosubflag)
  list(APPEND ${PROJECT_NAME}_flags $<$<COMPILE_LANGUAGE:Fortran>:-Wno-do-subscript>)
endif()

# add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-ffpe-trap=invalid,zero,overflow>")#,underflow)

# NOTE: don't use -march=native as GCC doesn't support all CPU arches with that option.
# add_compile_options(-mtune=native)

add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fimplicit-none>)

# --- IMPORTANT
add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-Werror=array-bounds;-fcheck=all>")
# --- IMPORTANT: options help trap array indexing/bounds errors at runtime

add_compile_options(-Wall)
# Fortran -Wall causes false positives in many Fortran projects
add_compile_options(
$<$<COMPILE_LANGUAGE:Fortran>:-Wno-uninitialized>
$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<VERSION_LESS:$<Fortran_COMPILER_VERSION>,10>>:-Wno-conversion>
)

if(dev)
  add_compile_options(-Wextra)
  # -Wpedantic makes too many false positives
endif()

# avoid backtrace that's unusable without -g
add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-fno-backtrace>)

# Wdo-subscript is known to warn on obvious non-problems
check_fortran_compiler_flag(-Wdo-subscript dosubflag)
if(dosubflag)
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Wno-do-subscript>)
endif()

# add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-ffpe-trap=invalid,zero,overflow>")#,underflow)

# lot of spurious warnings on allocatable scalar character
add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_VERSION:9.3.0>>:-Wno-maybe-uninitialized>)

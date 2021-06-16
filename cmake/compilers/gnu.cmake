# NOTE: don't use -march=native as GCC doesn't support all CPU arches with that option.
add_compile_options(-mtune=native)

# keep Wall in cmake_fortran_flags so FetchContent packages can override
# and thereby avoid useless megabytes of other project warnings
string(APPEND CMAKE_Fortran_FLAGS " -fimplicit-none")
if(dev)
  string(APPEND CMAKE_Fortran_FLAGS " -Wall -Wextra")
  # -Wpedantic makes too many false positives
else(dev)
  string(APPEND CMAKE_Fortran_FLAGS " -Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function")
endif(dev)

# avoid backtrace that's unusable without -g
if(CMAKE_BUILD_TYPE STREQUAL Release)
  string(APPEND CMAKE_Fortran_FLAGS " -fno-backtrace")
endif()

# Wdo-subscript is known to warn on obvious non-problems
check_compiler_flag(Fortran -Wdo-subscript dosubflag)
if(dosubflag)
  string(APPEND CMAKE_Fortran_FLAGS " -Wno-do-subscript")
endif()

string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -Werror=array-bounds -fcheck=all")
  # string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -ffpe-trap=invalid,zero,overflow")#,underflow)

# enforce Fortran 2018 standard
# this flag is buggy at least through GCC 9, causing fake "Common block" warnings
# check_compiler_flag(Fortran -std=f2018 f18flag)
# if(f18flag)
#   string(APPEND CMAKE_Fortran_FLAGS " -std=f2018")
# endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 9.3.0)
  # makes a lot of spurious warnngs on alloctable scalar character
  string(APPEND CMAKE_Fortran_FLAGS " -Wno-maybe-uninitialized")
endif()

check_compiler_flag(Fortran -fallow-argument-mismatch allow_mismatch_args)
if(allow_mismatch_args)
  set(gcc10opts -fallow-argument-mismatch)
endif()

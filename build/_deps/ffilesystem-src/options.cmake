option(ffilesystem_cpp "Use C++ filesystem for full functionality" on)
option(ffilesystem_fortran "use the Fortran interaces to C functions" on)
option(ffilesystem_cli "Build CLI" ${ffilesystem_IS_TOP_LEVEL})
option(ffilesystem_fallback "Fallback to non-C++ if C++ stdlib is not working" on)
option(ffilesystem_trace "debug trace output" off)
option(ffilesystem_locale "enable locale-based functions")
option(ffilesystem_extra "enable extra functions not strictly filesystem-based" on)
option(ffilesystem_unicode "Windows Unicode support" on)
option(ffilesystem_pic "Build position-independent code")

option(BUILD_SHARED_LIBS "Build shared libraries")
option(ffilesystem_tidy "Run clang-tidy on the code")
option(ffilesystem_cppcheck "Run cppcheck on the code")
option(ffilesystem_iwyu "Run include-what-you-use on the code")

option(ffilesystem_iso_fortran_binding "Use ISO_Fortran_binding.h for Cpp-Fortran interfaces")

option(ffilesystem_ENABLE_RPATH "Enable RPATH handling for transitive dependencies" on)

option(ffilesystem_BUILD_TESTING "Build tests")

if(NOT DEFINED boost_ut_url)
  set(boost_ut_url "https://github.com/boost-ext/ut/archive/59a9beba0763dbb45b3cc68e4cf484c659319a97.tar.gz")
endif()

file(GENERATE OUTPUT .gitignore CONTENT "*")

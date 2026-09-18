# Check compiler C++ capabilities
include(CheckSourceCompiles)

include_guard()

function(cpp_check)

set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

# normally this check is not triggered, it's a safety belt and more detailed approaches
# are possible for old compilers.
# See: https://github.com/gulrak/filesystem/blob/b1982f06c84f08a99fb90bac43c2d03712efe921/include/ghc/filesystem.hpp#L215
if(CMAKE_CXX_STANDARD LESS 17)
  check_cxx_symbol_exists(__cpp_lib_string_view "string_view" cpp_string_view)
  if(NOT cpp_string_view)
    message(WARNING "C++ stdlib string_view not detected in libstdc++ ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD}")
  endif()
endif()

set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)

# some compilers claim to have filesystem, but their libstdc++ doesn't have it.
set(_src "#include <filesystem>
int main () {
std::filesystem::path tgt;
return tgt.empty() ? 0 : 1;
}")
check_source_compiles(CXX "${_src}" HAVE_CXX_FILESYSTEM)

if(NOT HAVE_CXX_FILESYSTEM)
  message(STATUS "${PROJECT_NAME}: CMake ${CMAKE_VERSION} ${CMAKE_GENERATOR} proceeding with C++${CMAKE_CXX_STANDARD} fallback to <filesystem>: libstdc++ ${ffilesystem_stdcpp_version}  compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}  linker: ${CMAKE_CXX_COMPILER_LINKER_ID} ${CMAKE_CXX_COMPILER_LINKER_VERSION}
${HAVE_CXX_FILESYSTEM}
  see CMake log file ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeConfigureLog.yaml for details on variable HAVE_CXX_FILESYSTEM.")
  return()
endif()

if(CMAKE_CXX_STANDARD GREATER_EQUAL 20)

# for fs_get_modtime
check_source_compiles(CXX
"#include <chrono>
#include <filesystem>

int main(){

std::filesystem::file_time_type t_fs;
auto t_sys = std::chrono::clock_cast<std::chrono::system_clock>(t_fs);
return 0;
}"
ffilesystem_HAVE_CLOCK_CAST
)
else()
  unset(ffilesystem_HAVE_CLOCK_CAST CACHE)
endif()


if(ffilesystem_trace)
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

  check_cxx_symbol_exists(__cpp_lib_format "" cpp20_format)
  check_cxx_symbol_exists(__cpp_lib_ranges "" cpp20_ranges)

  check_cxx_symbol_exists(__cpp_lib_starts_ends_with "string" cpp20_string_ends_with)
  check_cxx_symbol_exists(__cpp_using_enum "" cpp20_using_enum)
  check_cxx_symbol_exists(__cpp_deduction_guides "" cpp17_deduction_guides)

  check_cxx_symbol_exists(__has_cpp_attribute "" cpp20_has_cpp_attribute)

  if(cpp20_has_cpp_attribute)

  check_source_compiles(CXX
  "#if !__has_cpp_attribute(likely)
  #error \"no likely attribute\"
  #endif"
  cpp20_likely
  )

  endif()
endif()

endfunction()

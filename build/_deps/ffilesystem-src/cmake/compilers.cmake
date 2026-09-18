include(CheckIncludeFile)
include(CheckSymbolExists)
include(CheckCXXSymbolExists)
include(CheckSourceCompiles)

# --- some compilers require these manual settings
unset(CMAKE_REQUIRED_FLAGS)
unset(CMAKE_REQUIRED_LIBRARIES)
unset(CMAKE_REQUIRED_DEFINITIONS)

# also MinGW Flang on ARM
#if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24 AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_GENERATOR STREQUAL "Unix Makefiles")
# otherwise failed to link since -lc++ is missing
set(ffilesystem_linker_lang)
if(LINUX AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM|NVHPC")
  # IntelLLVM|NVHPC need Fortran. For other compilers best to leave default.
  set(ffilesystem_linker_lang Fortran)
endif()

if(ffilesystem_cpp)
  include(${CMAKE_CURRENT_LIST_DIR}/CppCheck.cmake)
  cpp_check()
else()
  unset(HAVE_CXX_FILESYSTEM CACHE)
endif()

if(UNIX AND BUILD_SHARED_LIBS)
  block()
  list(APPEND CMAKE_REQUIRED_LIBRARIES ${CMAKE_DL_LIBS})

  set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
  check_cxx_symbol_exists(dladdr "dlfcn.h" ffilesystem_HAVE_DLADDR)

  endblock()
endif()

# --- deeper filesystem check: C, C++ and Fortran compiler ABI compatibility

if(HAVE_CXX_FILESYSTEM)

if(ffilesystem_trace)
  # this is not on by default because it gave a false negative on some Cray systems
  # while bypassing it let all of Ffilesystem work.
  include(${CMAKE_CURRENT_LIST_DIR}/FScheck.cmake)
endif()

elseif(LINUX OR ANDROID)
  block()
  set(CMAKE_REQUIRED_DEFINITIONS -D_DEFAULT_SOURCE)

  check_symbol_exists(copy_file_range "sys/types.h;unistd.h" ffilesystem_HAVE_COPY_FILE_RANGE)

  check_symbol_exists(statx "sys/stat.h;fcntl.h" ffilesystem_HAVE_STATX)
  endblock()
endif()


if(WIN32)
  check_cxx_symbol_exists(GetFileInformationByName "Windows.h" ffilesystem_HAVE_GetFileInformationByName)
elseif(LINUX OR ANDROID)
  block()
  set(CMAKE_REQUIRED_DEFINITIONS -D_DEFAULT_SOURCE)
  check_cxx_symbol_exists(statx "sys/stat.h;fcntl.h" ffilesystem_HAVE_STATX)
  endblock()
endif()


if(ffilesystem_cpp AND NOT ffilesystem_fallback AND NOT HAVE_CXX_FILESYSTEM)
  message(FATAL_ERROR "C++ filesystem not available. To fallback to plain C++:
  cmake -Dffilesystem_fallback=on -B build"
  )
endif()

# --- END COMPILER CHECKS


# --- C / C++ compile flags
if(NOT MSVC AND CMAKE_C_COMPILER_ID MATCHES "Clang|GNU|^Intel")
  target_compile_options(ffilesystem PRIVATE
  "$<$<AND:$<COMPILE_LANGUAGE:C,CXX>,$<CONFIG:Debug,RelWithDebInfo>>:-Wall;-Wextra>"
  "$<$<COMPILE_LANGUAGE:C>:-Werror=implicit-function-declaration>"
  )
endif()

if(MSVC OR CMAKE_CXX_SIMULATE_ID STREQUAL "MSVC")
  # the CMAKE_CXX_SIMULATE_ID catches when frontend is GNU but MSVC runtime headers are used
  target_compile_definitions(ffilesystem PRIVATE "$<$<COMPILE_LANGUAGE:C,CXX>:_CRT_SECURE_NO_WARNINGS>")

  # if, not elseif, because IntelLLVM uses MSVC flags
  # /wd4996 quiets warning: 'GetVersionExA' is deprecated
  target_compile_options(ffilesystem PRIVATE
  "$<$<COMPILE_LANGUAGE:C,CXX>:/W3;/wd4996>"
  "$<$<AND:$<COMPILE_LANGUAGE:C,CXX>,$<CXX_COMPILER_ID:Clang>>:-Wno-deprecated-declarations>"
  )
  if(ffilesystem_unicode)
    add_compile_definitions("$<$<COMPILE_LANGUAGE:C,CXX>:_UNICODE>")
  endif()
endif()

target_compile_options(ffilesystem PRIVATE "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang,GNU>:-Wold-style-cast>")

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.14)
  # MSVC has __cpluscplus = 199711L by default, which is C++98!
  # MSVC 15.7 (19.14) added /Zc:__cplusplus to set __cplusplus to the true value.
  # oneAPI since 2023.1 sets __cplusplus to the true value with MSVC by auto-setting this flag.
  target_compile_options(ffilesystem PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:/Zc:__cplusplus>")
endif()

# too many false positives esp. regarding argv[1] and higher
# target_compile_options(ffilesystem PRIVATE "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang,IntelLLVM>:-Wunsafe-buffer-usage>")


if(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
  target_compile_options(ffilesystem PRIVATE "$<$<AND:$<COMPILE_LANGUAGE:C,CXX>,$<CONFIG:Debug>>:-Rno-debug-disables-optimization>")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10")
  target_compile_options(ffilesystem PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wno-attributes>")
  # this is for UNLIKELY/LIKELY macros
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
  target_compile_options(ffilesystem PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:--diag_suppress=code_is_unreachable>")
endif()

# --- Fortran compile flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")

target_compile_options(ffilesystem PRIVATE
"$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-warn>"
"$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-traceback;-check;-debug>"
)

# this flag needs to be applied EVERYWHERE incl. submodule projects with add_compile_options()
# or runtime errors / weird behavior with unresolved procedures that actually exist.
# -standard-semantics is no good because it breaks linkage within oneAPI itself e.g. oneMPI library!
if(WIN32)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:/fpscomp:logicals>")
else()
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-fpscomp;logicals>")
endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

target_compile_options(ffilesystem PRIVATE
"$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-Wall;-Wextra>"
"$<$<COMPILE_LANGUAGE:Fortran>:-fimplicit-none>"
"$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-fcheck=all;-Werror=array-bounds>"
"$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-fno-backtrace>"
)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  # C_BOOL correctness
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Munixlogical>")
endif()

# --- code coverage
if(ffilesystem_coverage)
  if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?[Cc]lang")
      message(STATUS "${PROJECT_NAME} Clang coverage")
      target_compile_options(ffilesystem PRIVATE -fprofile-instr-generate -fcoverage-mapping)
      target_link_options(ffilesystem PRIVATE -fprofile-instr-generate)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      message(STATUS "${PROJECT_NAME} GCC coverage")
      target_compile_options(ffilesystem PRIVATE --coverage -O0 -g)
      target_link_options(ffilesystem PRIVATE --coverage)
  else()
      message(FATAL_ERROR "Code coverage is only supported with Clang or GCC")
  endif()

  # Optional: custom target that runs tests + coverage collection
  add_custom_target(coverage
      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --test-dir ${PROJECT_BINARY_DIR}
      COMMAND ${CMAKE_COMMAND} -E echo "Generating coverage report..."
      COMMAND ${CMAKE_COMMAND} -E env GCOV_PREFIX=${PROJECT_BINARY_DIR}/coverage gcovr -r ${CMAKE_SOURCE_DIR} --html --html-details -o ${PROJECT_BINARY_DIR}/coverage/coverage.html
  )
endif()

# --- clang-tidy
if(ffilesystem_tidy)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy REQUIRED
  PATHS /opt/homebrew/opt/llvm/bin
  )
  set(tidy_cmd ${CLANG_TIDY_EXE} -format-style=file)
  set(CMAKE_C_CLANG_TIDY ${tidy_cmd})
  set(CMAKE_CXX_CLANG_TIDY ${tidy_cmd})
endif()

# --- IWYU
if(ffilesystem_iwyu)
  find_program(IWYU_EXE NAMES include-what-you-use REQUIRED)
  message(STATUS "IWYU_EXE: ${IWYU_EXE}")
  set(iwyu_cmd ${IWYU_EXE})
  set(CMAKE_C_INCLUDE_WHAT_YOU_USE ${iwyu_cmd})
  set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE ${iwyu_cmd})
endif()

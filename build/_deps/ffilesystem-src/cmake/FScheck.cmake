block()

set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)

if(NOT DEFINED ffilesystem_abi_ok)

message(CHECK_START "checking that compilers can link C++ Filesystem together")

try_compile(ffilesystem_abi_ok
PROJECT fs_check
SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/fs_check
CMAKE_FLAGS
  -DGNU_stdfs=${GNU_stdfs}
  -Dffilesystem_linker_lang=${ffilesystem_linker_lang}
  -Dffilesystem_fortran:BOOL=${ffilesystem_fortran}
)

if(ffilesystem_abi_ok)
  message(CHECK_PASS "OK")
else()
  message(CHECK_FAIL "Failed")
  # message(WARNING "Disabling C++ filesystem due to ABI-incompatible compilers or compiler problem")
  # set(HAVE_CXX_FILESYSTEM false CACHE BOOL "ABI problem with C++ filesystem" FORCE)
endif()

endif()

endblock()

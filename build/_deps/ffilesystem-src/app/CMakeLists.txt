if(CMAKE_CXX_STANDARD LESS 20 AND "cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  set(CMAKE_CXX_STANDARD 20)
endif()

add_executable(ffilesystem_bench benchmark.cpp)
target_link_libraries(ffilesystem_bench PRIVATE ffilesystem)

add_executable(fs_cli main.cpp $<$<BOOL:${WIN32}>:windows_long_paths.manifest>)
target_link_libraries(fs_cli PRIVATE ffilesystem)
target_compile_definitions(fs_cli PRIVATE
"$<$<BOOL:${HAVE_CXX_FILESYSTEM}>:HAVE_CXX_FILESYSTEM>"
"$<$<BOOL:${ffilesystem_extra}>:ffilesystem_extra>"
)
install(TARGETS fs_cli EXPORT ffilesystem-targets
BUNDLE DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(libc_version libc.cpp)
add_executable(libcpp_version libcxx.cpp)

if(ffilesystem_fortran)
  add_executable(filesystem_cli fortran/main.f90 $<$<BOOL:${WIN32}>:windows_long_paths.manifest>)
  target_link_libraries(filesystem_cli PRIVATE ffilesystem)
  set_property(TARGET filesystem_cli PROPERTY LINKER_LANGUAGE ${ffilesystem_linker_lang})

  install(TARGETS filesystem_cli EXPORT ffilesystem-targets)
endif()

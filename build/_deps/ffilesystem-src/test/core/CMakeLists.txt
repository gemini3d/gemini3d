set_property(DIRECTORY PROPERTY LABELS "core")

if(ffilesystem_fortran)

set(tests is_absolute absolute
canonical realpath
env
equivalent
exe
expanduser
filename
filesize space
is_dir
is_empty
is_subdir join
memory
normal parent root
resolve
relative proximate
reserved safe
set_cwd
separator
stem
suffix
with_suffix
owner
permissions
copyfile
mkdir
touch
which
)

if(HAVE_F03TYPE)
  list(APPEND tests path_t)
endif()

if(CYGWIN)
  list(APPEND tests cygwin windows)
elseif(WIN32)
  list(APPEND tests windows)
endif()

foreach(t IN LISTS tests)

add_executable(test_${t}_fortran test_${t}.f90)
target_link_libraries(test_${t}_fortran PRIVATE ffilesystem)
target_compile_options(test_${t}_fortran PRIVATE ${ffilesystem_fortran_test_flags})
set_property(TARGET test_${t}_fortran PROPERTY LINKER_LANGUAGE ${ffilesystem_linker_lang})

add_test(NAME Fortran_${t} COMMAND test_${t}_fortran)
set_tests_properties(Fortran_${t} PROPERTIES
LABELS Fortran
SKIP_RETURN_CODE 77
WORKING_DIRECTORY $<TARGET_FILE_DIR:test_${t}_fortran>
)

endforeach()

set_tests_properties(Fortran_relative Fortran_proximate PROPERTIES DISABLED $<NOT:$<BOOL:${HAVE_CXX_FILESYSTEM}>>)

if(ffilesystem_HAVE_ISO_FORTRAN_BINDING)
  add_executable(test_realpath_cfi_fortran test_realpath_cfi.f90)
  target_link_libraries(test_realpath_cfi_fortran PRIVATE ffilesystem)
  target_compile_options(test_realpath_cfi_fortran PRIVATE ${ffilesystem_fortran_test_flags})
  set_property(TARGET test_realpath_cfi_fortran PROPERTY LINKER_LANGUAGE ${ffilesystem_linker_lang})

  add_test(NAME Fortran_realpath_cfi COMMAND test_realpath_cfi_fortran)
  set_tests_properties(Fortran_realpath_cfi PROPERTIES
    LABELS "Fortran;CFI"
    WORKING_DIRECTORY $<TARGET_FILE_DIR:test_realpath_cfi_fortran>
  )
else()
  add_custom_target(test_realpath_cfi_fortran
    COMMAND ${CMAKE_COMMAND} -E echo "Skipping test_realpath_cfi_fortran: ISO_Fortran_binding.h not found or working"
  )
endif()

endif()


# --- Boost UT tests
if(HAVE_BOOST_UT)

set(boost_ut_tests
absolute
canonical
copyfile
drop_slash
empty
env
equivalent
exe
expanduser
filename
filesize
hard_link
is_char
is_empty
is_prefix
is_subdir
join
memory
mkdir
normal
OnDisk
owner
parent
permissions
pure
random
remove
rename
reserved
root
safe
separator
set_cwd
space
system
utf8
which
with_suffix
)

if(NOT ANDROID)
  list(APPEND boost_ut_tests is_fifo)
endif()

if(WIN32 OR CYGWIN)
  list(APPEND boost_ut_tests is_appexec_alias windows)
endif()

if(HAVE_CXX_FILESYSTEM)
  list(APPEND boost_ut_tests proximate relative)
endif()
foreach(t IN LISTS boost_ut_tests)
  add_executable(test_${t}_cpp test_${t}.cpp)
  target_link_libraries(test_${t}_cpp PRIVATE ffilesystem)
  target_include_directories(test_${t}_cpp PRIVATE ${boost_ut_SOURCE_DIR}/include ${ffilesystem_SOURCE_DIR}/test)
  target_compile_options(test_${t}_cpp PRIVATE ${ffilesystem_cpp_test_flags})

  add_test(NAME ${t}_cpp COMMAND test_${t}_cpp)
  set_tests_properties(${t}_cpp PROPERTIES
  LABELS "Cpp"
  WORKING_DIRECTORY $<TARGET_FILE_DIR:test_${t}_cpp>
  )
endforeach()

endif()

# C

foreach(t IN ITEMS overflow)

add_executable(test_${t}_c test_${t}.c)
target_link_libraries(test_${t}_c PRIVATE ffilesystem)

add_test(NAME C_${t} COMMAND test_${t}_c)
set_tests_properties(C_${t} PROPERTIES
LABELS "C"
WORKING_DIRECTORY $<TARGET_FILE_DIR:test_${t}_c>
)

endforeach()

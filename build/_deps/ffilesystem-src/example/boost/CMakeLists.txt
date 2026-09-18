cmake_minimum_required(VERSION 3.30...4.1)

project(BoostFilesystem LANGUAGES CXX)

enable_testing()

find_package(Boost REQUIRED COMPONENTS filesystem)

add_executable(windows_root boost_root.cpp)
target_link_libraries(windows_root PRIVATE Boost::filesystem)
target_compile_features(windows_root PRIVATE cxx_std_17)

add_test(NAME WindowsRootBoost COMMAND windows_root)

if(MSVC)
  set(boost_bindir "${Boost_INCLUDE_DIRS}/../../../installed/x64-windows/bin")
  if(IS_DIRECTORY ${boost_bindir})
    file(REAL_PATH "${boost_bindir}" boost_bindir)
    set_tests_properties(WindowsRootBoost PROPERTIES ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${boost_bindir}")
  endif()
endif()

file(GENERATE OUTPUT .gitignore CONTENT "*")

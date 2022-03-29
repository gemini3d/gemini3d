# package for stubs
function(stub_package name)

include(CMakePackageConfigHelpers)

configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/package/${name}-config.cmake.in
${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${name}-config.cmake
INSTALL_DESTINATION cmake
)

write_basic_package_version_file(
${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${name}-config-version.cmake
COMPATIBILITY SameMinorVersion
)

install(EXPORT ${name}-targets
NAMESPACE ${name}::
DESTINATION cmake
)

install(FILES
${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${name}-config.cmake
${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${name}-config-version.cmake
DESTINATION cmake
)

endfunction(stub_package)

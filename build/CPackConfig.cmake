# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_BUILD_SOURCE_DIRS "/home/runner/work/gemini3d/gemini3d;/home/runner/work/gemini3d/gemini3d/build")
set(CPACK_CMAKE_GENERATOR "Ninja")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libhdf5-dev (>=1.10)")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/tmp/cmake325/cmake-3.25.3-linux-x86_64/share/cmake-3.25/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "gemini3d built using CMake")
set(CPACK_GENERATOR "TBZ2")
set(CPACK_INSTALL_CMAKE_PROJECTS "/home/runner/work/gemini3d/gemini3d/build;gemini3d;ALL;/")
set(CPACK_INSTALL_PREFIX "/home/runner/work/gemini3d/gemini3d/build/local")
set(CPACK_MODULE_PATH "/home/runner/work/gemini3d/gemini3d/build/_deps/h5fortran-src/cmake")
set(CPACK_NSIS_DISPLAY_NAME "gemini3d 2.0.0")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_PACKAGE_NAME "gemini3d 2.0.0")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_OBJCOPY_EXECUTABLE "/usr/bin/objcopy")
set(CPACK_OBJDUMP_EXECUTABLE "/usr/bin/objdump")
set(CPACK_OUTPUT_CONFIG_FILE "/home/runner/work/gemini3d/gemini3d/build/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/tmp/cmake325/cmake-3.25.3-linux-x86_64/share/cmake-3.25/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "3-D ionospheric model")
set(CPACK_PACKAGE_FILE_NAME "gemini3d-2.0.0-Linux")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/gemini3d/gemini")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "gemini3d 2.0.0")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "gemini3d 2.0.0")
set(CPACK_PACKAGE_NAME "gemini3d")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "SciVision")
set(CPACK_PACKAGE_VERSION "2.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_READELF_EXECUTABLE "/usr/bin/readelf")
set(CPACK_RESOURCE_FILE_LICENSE "/home/runner/work/gemini3d/gemini3d/build/_deps/h5fortran-src/LICENSE")
set(CPACK_RESOURCE_FILE_README "/home/runner/work/gemini3d/gemini3d/build/_deps/h5fortran-src/README.md")
set(CPACK_RESOURCE_FILE_WELCOME "/tmp/cmake325/cmake-3.25.3-linux-x86_64/share/cmake-3.25/Templates/CPack.GenericWelcome.txt")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_IGNORE_FILES ".git/;.github/;.vscode/;.mypy_cache/;_CPack_Packages/;/home/runner/work/gemini3d/gemini3d/build/;/home/runner/work/gemini3d/gemini3d/build/_deps/h5fortran-build/;.archive/;concepts/;paper/")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/runner/work/gemini3d/gemini3d/build/CPackSourceConfig.cmake")
set(CPACK_SYSTEM_NAME "Linux")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "Linux")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/runner/work/gemini3d/gemini3d/build/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()

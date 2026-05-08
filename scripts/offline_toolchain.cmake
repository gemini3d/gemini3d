# use this to build Gemini3D using offline libraries like:
#
#    cmake -P scripts/offline_libraries.cmake
#
#    cmake --toolchain scripts/offline_toolchain.cmake -B build
#
#    cmake --build build
#
cmake_minimum_required(VERSION 3.21)

set(offline_dir ${CMAKE_CURRENT_LIST_DIR}/../offline)

set(FETCHCONTENT_SOURCE_DIR_FFILESYSTEM ${offline_dir}/ffilesystem-src)
set(FETCHCONTENT_SOURCE_DIR_GLOW ${offline_dir}/glow-src)
set(FETCHCONTENT_SOURCE_DIR_H5FORTRAN ${offline_dir}/h5fortran-src)
set(FETCHCONTENT_SOURCE_DIR_HWM14 ${offline_dir}/hwm14-src)
set(FETCHCONTENT_SOURCE_DIR_MSIS ${offline_dir}/msis-src)

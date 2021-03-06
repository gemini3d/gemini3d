cmake_minimum_required(VERSION 3.20...3.21)

set(version 1.5.0)

set(prefix "~/zstd-${version}")

set(CMAKE_TLS_VERIFY true)

if(DEFINED ENV{TMPDIR})
  set(tmpdir $ENV{TMPDIR})
elseif(IS_DIRECTORY /var/tmp)
  set(tmpdir /var/tmp)
elseif(IS_DIRECTORY /tmp)
  set(tmpdir /tmp)
else()
  set(tmpdir ~/tmp)
endif()

if(CMAKE_VERSION VERSION_LESS 3.21)
  get_filename_component(tmpdir ${tmpdir} ABSOLUTE)
else()
  file(REAL_PATH ${tmpdir} tmpdir EXPAND_TILDE)
endif()

set(name zstd-${version}.tar.gz)
set(archive ${tmpdir}/${name})

set(src ${tmpdir}/zstd-${version}/build/cmake)
set(build ${src}/build)

if(NOT IS_DIRECTORY ${src})
  file(DOWNLOAD https://github.com/facebook/zstd/releases/download/v${version}/${name} ${archive}
    INACTIVITY_TIMEOUT 15)
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${tmpdir})
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX:PATH=${prefix} -S ${src} -B ${build}
COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND ${CMAKE_COMMAND} --build ${build} --parallel
COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND ${CMAKE_COMMAND} --install ${build}
COMMAND_ERROR_IS_FATAL ANY)

message(STATUS "Please add ${CMAKE_INSTALL_PREFIX}/bin to environment variable PATH")

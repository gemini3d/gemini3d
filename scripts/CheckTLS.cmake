cmake_minimum_required(VERSION 3.20)

function(check_tls)
# some CMake may not have SSL/TLS enabled, or may have missing/broken system certificates.
# this is a publicly-usable service (as per their TOS)

set(url https://www.howsmyssl.com/a/check)
cmake_path(APPEND tls_check_file ${PROJECT_BINARY_DIR} check_tls.json)

if(EXISTS ${tls_check_file})
  # we've already checked this computer's TLS configuration.
  return()
endif()

message(STATUS "CheckTLS: ${url} => ${tls_check_file}")
file(DOWNLOAD ${url} ${tls_check_file} INACTIVITY_TIMEOUT 10)

file(READ ${tls_check_file} json)

string(JSON rating ERROR_VARIABLE e GET ${json} rating)

message(STATUS "TLS status: ${rating}")
if(NOT rating STREQUAL "Probably Okay")
  message(WARNING "TLS seems to be broken on your system. Download will probably fail.  ${rating}")
endif()

endfunction(check_tls)

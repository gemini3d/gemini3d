# use GNU Patch from any platform
#
# Functions
# ---------
#
# patch_file(in_file patch_file)
#     Apply patch_file to in_file via GNU Patch.

if(WIN32)
  find_package(Msys)
  if(MSYS_INSTALL_PATH)
    find_program(PATCH
    NAMES patch
    HINTS ${MSYS_INSTALL_PATH}
    PATH_SUFFIXES bin usr/bin
    )
  endif()

  if(NOT PATCH)
    find_program(WSL NAMES wsl)
  endif()
else()
  find_program(PATCH NAMES patch)
endif()


function(patch_file in_file patch_file)

if(PATCH)
  execute_process(COMMAND ${PATCH} ${in_file} ${patch_file}
  TIMEOUT 15
  COMMAND_ERROR_IS_FATAL ANY
  )
elseif(WSL)
  execute_process(COMMAND ${WSL} wslpath ${in_file}
  TIMEOUT 5
  OUTPUT_VARIABLE in_wsl
  COMMAND_ERROR_IS_FATAL ANY
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(COMMAND ${WSL} wslpath ${patch_file}
  TIMEOUT 5
  OUTPUT_VARIABLE patch_wsl
  COMMAND_ERROR_IS_FATAL ANY
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(COMMAND ${WSL} patch ${in_wsl} ${patch_wsl}
  TIMEOUT 15
  COMMAND_ERROR_IS_FATAL ANY
  )
else()
  message(FATAL_ERROR "Could not find patch program")
endif()
endfunction(patch_file)

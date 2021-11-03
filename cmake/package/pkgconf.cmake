# --- generate pkg-config .pc

set(pc_requires "h5fortran lapack")

set(pc_filename ${PROJECT_NAME}.pc)
configure_file(${CMAKE_CURRENT_LIST_DIR}/pkgconf.pc.in ${pc_filename} @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${pc_filename} DESTINATION pkgconfig)

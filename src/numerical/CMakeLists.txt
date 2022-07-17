# --- generate constants
configure_file(constants/phys_consts.in.f90 phys_consts.f90 @ONLY)
add_library(const OBJECT ${CMAKE_CURRENT_BINARY_DIR}/phys_consts.f90)

# --- interpolation
add_subdirectory(interpolation)

add_subdirectory(grid)

add_subdirectory(calculus)

add_subdirectory(diffusion)

add_subdirectory(advection)

add_subdirectory(potential)

add_subdirectory(coord)

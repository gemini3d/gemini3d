set_directory_properties(PROPERTIES LABELS interpolation)

add_library(interp OBJECT interpolation.f90 interp2d.f90)
target_link_libraries(interp PRIVATE const)

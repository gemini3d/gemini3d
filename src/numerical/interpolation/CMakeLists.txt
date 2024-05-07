set_property(DIRECTORY PROPERTY  LABELS interpolation)

add_library(interp OBJECT interpolation.f90 interp2d.f90)
target_link_libraries(interp PRIVATE const)

# apply compiler options
get_property(tgts DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
foreach(t IN LISTS tgts)
  target_compile_options(${t} PRIVATE ${gemini3d_flags})
endforeach()

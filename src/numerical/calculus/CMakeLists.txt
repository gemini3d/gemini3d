configure_file(calculus.in.f90 calculus.f90 @ONLY)
add_library(calculus OBJECT
${CMAKE_CURRENT_BINARY_DIR}/calculus.f90 gradient.f90 integral.f90 div.f90
)
target_link_libraries(calculus PRIVATE meshobj)

# apply compiler options
get_property(tgts DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
foreach(t IN LISTS tgts)
  target_compile_options(${t} PRIVATE ${gemini3d_flags})
endforeach()

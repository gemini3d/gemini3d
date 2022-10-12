set_property(DIRECTORY PROPERTY LABELS diffusion)

set(_tmpdifffn ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion1d.h5)

add_executable(test_diffusion1d test_diffusion1D.f90)
# lapack needed for MacOS
target_link_libraries(test_diffusion1d PRIVATE PDEparabolic gbsv const
LAPACK::LAPACK
h5fortran::h5fortran
)

# WORKING_DIRECTORY is needed to write output file to repeatable path.
add_test(NAME diffusion1
COMMAND test_diffusion1d ${_tmpdifffn}
WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)
set_tests_properties(diffusion1 PROPERTIES
LABELS unit
FIXTURES_SETUP gemini_diffusion_fxt
)
dll_test_path("h5fortran::h5fortran;ZLIB::ZLIB;LAPACK::LAPACK" diffusion1)

add_test(NAME diffusion1:python
COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_diffusion1D.py ${_tmpdifffn}
)
set_tests_properties(diffusion1:python PROPERTIES
LABELS "unit;python"
REQUIRED_FILES ${_tmpdifffn}
FIXTURES_REQUIRED gemini_diffusion_fxt
DISABLED $<NOT:$<BOOL:${H5PY_FOUND}>>
)

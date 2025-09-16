set_property(DIRECTORY PROPERTY LABELS potential)

# --- potential 2D test

set(_potential2d_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential2d.h5)

add_executable(gemini_potential2d test_potential2d.f90 $<TARGET_OBJECTS:PDEelliptic>)
set_target_properties(gemini_potential2d PROPERTIES LABELS "unit")
target_link_libraries(gemini_potential2d PRIVATE const
h5fortran::h5fortran MUMPS::MUMPS MPI::MPI_Fortran
)
target_include_directories(gemini_potential2d PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
# this include is for robustness so mpi_f08.mod isn't missed e.g. on Linux systems

# LINK_INTERFACE_MULTIPLICITY and IMPORTED_LINK_INTERFACE_MULTIPLICITY do not help here.
# adding INTERFACE to Mumps project didn't help either
# nor did adding them in mumps.cmake
# it seems to have to be right here.
# gemini_potential2d has been an intermittent source of link-order issues,
# where mpi gets linked in the middle instead of at the end.
# the TARGET_OBJECTS etc in this directory are a careful effort to avoid this problem

add_test(NAME GeminiPotential2d COMMAND gemini_potential2d ${_potential2d_testfile})
test_mpi_launcher(gemini_potential2d GeminiPotential2d ${Ncpu})

add_test(NAME potential2d:python
COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_potential2d.py ${_potential2d_testfile}
)
set_tests_properties(potential2d:python PROPERTIES
REQUIRED_FILES ${_potential2d_testfile}
LABELS "unit;python"
FIXTURES_REQUIRED gemini_potential_fxt
DISABLED $<NOT:$<BOOL:${H5PY_FOUND}>>
)


# --- potential 3D test

set(_potential3d_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential3d.h5)

add_executable(gemini_potential3d test_potential3d.f90 $<TARGET_OBJECTS:PDEelliptic>)
set_target_properties(gemini_potential3d PROPERTIES LABELS "unit")
target_link_libraries(gemini_potential3d PRIVATE const
h5fortran::h5fortran MUMPS::MUMPS MPI::MPI_Fortran
)
target_include_directories(gemini_potential3d PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
# this include is for robustness so mpi_f08.mod isn't missed e.g. on Linux systems

add_test(NAME GeminiPotential3d COMMAND gemini_potential3d ${_potential3d_testfile})
test_mpi_launcher(gemini_potential3d GeminiPotential3d ${Ncpu})

set_tests_properties(GeminiPotential2d GeminiPotential3d PROPERTIES
FIXTURES_SETUP gemini_potential_fxt
FIXTURES_REQUIRED mpi_fxt
RESOURCE_LOCK cpu_mpi
LABELS "unit"
)

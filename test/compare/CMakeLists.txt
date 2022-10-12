add_executable(gemini3d.compare)
target_sources(gemini3d.compare PRIVATE
gemini_compare.f90 compare.f90 compare_grid_h5.f90 compare_out_h5.f90 compare_in_h5.f90
)
target_link_libraries(gemini3d.compare PRIVATE
assert gemini3d_config reader timeutils const gemini3d_sysinfo
ffilesystem::filesystem
h5fortran::h5fortran
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_property(TARGET gemini3d.compare PROPERTY LINKER_LANGUAGE Fortran)
else()
  set_property(TARGET gemini3d.compare PROPERTY LINKER_LANGUAGE CXX)
endif()

set_property(TARGET gemini3d.compare PROPERTY RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set_property(TARGET gemini3d.compare PROPERTY DEBUG_POSTFIX .debug)
set_property(TARGET gemini3d.compare PROPERTY RELWITHDEBINFO_POSTFIX .debug)

# --- Install

install(TARGETS gemini3d.compare
EXPORT ${PROJECT_NAME}-targets
)


# --- helper functions

function(matlab_compare outdir refdir name)

add_test(NAME gemini:compare:${name}:matlab
COMMAND ${Matlab_MAIN_PROGRAM} -batch "gemini3d.compare('${outdir}', '${refdir}')"
WORKING_DIRECTORY ${matgemini_SOURCE_DIR}
)

set_tests_properties(gemini:compare:${name}:matlab PROPERTIES
FIXTURES_REQUIRED "${name}:run_fxt"
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
ENVIRONMENT "MATLABPATH=${MATLABPATH}"
LABELS "compare;matlab"
)

endfunction(matlab_compare)


function(python_compare outdir refdir name)

add_test(NAME gemini:compare:${name}:python
COMMAND ${Python_EXECUTABLE} -m gemini3d.compare ${outdir} ${refdir} -file_format h5)

set_tests_properties(gemini:compare:${name}:python PROPERTIES
FIXTURES_REQUIRED ${name}:run_fxt
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS "compare;python"
DISABLED $<NOT:$<BOOL:${PYGEMINI_DIR}>>
)

endfunction(python_compare)


function(fortran_compare outdir refdir name)

add_test(NAME gemini:compare:${name}
COMMAND gemini3d.compare ${outdir} ${refdir})

set_tests_properties(gemini:compare:${name} PROPERTIES
FIXTURES_REQUIRED ${name}:run_fxt
RESOURCE_LOCK $<$<BOOL:${WIN32}>:cpu_mpi>
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS compare
DISABLED $<NOT:$<TARGET_EXISTS:gemini3d.compare>>
)
dll_test_path("ffilesystem::filesystem;h5fortran::h5fortran;HDF5::HDF5" gemini:compare:${name})

# resource_lock compare for Windows, which can take 100x longer when run
# at same time with non-dependent sim runs.
# it's not a problem to run multiple compare at once, but it is a problem
# to run gemini3d.compare at same time as gemini.bin, even on different sims

endfunction(fortran_compare)


function(compare_gemini_output name outdir refdir)

if(matlab)
  matlab_compare(${outdir} ${refdir} ${name})
endif()

if(python)
  python_compare(${outdir} ${refdir} ${name})
endif()

fortran_compare(${outdir} ${refdir} ${name})

endfunction(compare_gemini_output)

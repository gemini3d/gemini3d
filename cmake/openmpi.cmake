# https://github.com/open-mpi/ompi/issues/7393
# https://github.com/gerlero/openfoam-app/pull/112
# https://apple.stackexchange.com/a/287710
# *** OpenMPI workaround
if(UNIX AND MPI_C_LIBRARY_VERSION_STRING MATCHES "Open[ ]?MPI")
  if(NOT DEFINED mpi_tmpdir)
    execute_process(COMMAND mktemp -d /tmp/mpi-XXXXXXXX
    RESULT_VARIABLE ret
    OUTPUT_VARIABLE mpi_tmpdir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT ret EQUAL 0)
      message(FATAL_ERROR "could not create MPI working dir via mktemp -d: ${ret}")
    endif()
    set(mpi_tmpdir ${mpi_tmpdir} CACHE PATH "MPI working dir")
    message(STATUS "${ret}: Created MPI working dir ${mpi_tmpdir}")
  endif()
endif()
# *** END OpenMPI workaround

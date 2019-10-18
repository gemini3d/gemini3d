test_mpi_exe = executable('test_mpi', 'test_mpi.f90',
  dependencies: mpi,
  link_language: 'fortran')
test('MPI_OK', test_mpi_exe,
  timeout: 10,
  suite: 'unit')


test_mumps_exe = executable('test_mumps', 'test_mumps.F90',
  dependencies: [mumps, mpi],
  fortran_args: real_bits)
test('MUMPS', mpiexec,
  args: ['-np', '2', test_mumps_exe],
  workdir: meson.current_source_dir(),
  timeout: 15,
  is_parallel: false,
  suite: 'unit')

test_mpi_exe = executable('test_mpi', 'test_mpi.f90',
  dependencies: mpi)
test('MPI_OK', test_mpi_exe,
  timeout: 5,
  suite: 'unit')


test_mumps_exe = executable('test_mumps', 'test_mumps.F90',
  dependencies: [mumps, mpi],
  fortran_args: real_bits)
test('MUMPS', mpiexec,
  args: ['-np', '1', test_mumps_exe],
  workdir: meson.current_source_dir(),
  timeout: 5,
  suite: 'unit')
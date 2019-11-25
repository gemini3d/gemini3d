test_mpi_exe = executable('test_mpi', 'test_mpi.f90',
  dependencies: mpi,
  link_language: 'fortran')
test('MPI basic', test_mpi_exe,
  timeout: 10,
  suite: 'unit')

test_scalapack_exe = executable('test_scalapack', 'test_scalapack.f90',
  dependencies: [scalapack, blacs, mpi],
  link_language: 'fortran')
test('Scalapack basic', test_scalapack_exe,
  timeout: 10,
  suite: 'unit')



test_mumps_exe = executable('test_mumps', 'test_mumps.F90',
  dependencies: [mumps, scalapack, blacs, mpi],
  fortran_args: real_bits,
  link_language: 'fortran')
test('MUMPS basic', mpiexec,
  args: ['-np', '2', test_mumps_exe],
  workdir: meson.current_source_dir(),
  timeout: 15,
  is_parallel: false,
  suite: 'unit')

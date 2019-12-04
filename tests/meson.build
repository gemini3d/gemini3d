test_mpi_exe = executable('test_mpi', 'test_mpi.f90',
  dependencies: mpi,
  link_language: 'fortran')
test('mpi_basic', test_mpi_exe,
  timeout: 10,
  env: nomalloc_env,
  suite: 'unit')


test_scalapack_exe = executable('test_scalapack', 'test_scalapack.f90',
  dependencies: [scalapack, blacs, lapack, mpi],
  link_language: 'fortran')
test('scalapack_basic', test_scalapack_exe,
  timeout: 10,
  env: nomalloc_env,
  suite: 'unit')

test_mumps_exe = executable('test_mumps', 'test_mumps.F90',
  dependencies: [mumps, scalapack, blacs, lapack, mpi],
  fortran_args: real_bits,
  link_language: 'fortran')
test('mumps_basic', mpiexec,
  args: ['-np', '2', test_mumps_exe, files('input_simpletest_real.txt')],
  is_parallel: false,
  env: nomalloc_env,
  suite: 'unit')

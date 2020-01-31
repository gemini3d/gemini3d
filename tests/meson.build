test_mpi_exe = executable('test_mpi', 'test_mpi.f90',
  dependencies: mpi)
test('mpi_basic', test_mpi_exe,
  timeout: 10,
  env: nomalloc_env,
  suite: 'unit')


if arith.contains('d')
test_scalapack_d = executable('test_scalapack_d', 'test_scalapack_d.f90',
  dependencies: [scalapack, blacs, lapack, mpi])
test('scalapack_d', test_scalapack_d,
  timeout: 10,
  env: nomalloc_env,
  suite: 'unit')
endif

test_mumps_exe = executable('test_mumps',
  sources : ['test_mumps.f90', mumpsifc_src],
  dependencies: [mumps, scalapack, blacs, lapack, mpi])
test('mumps_basic', mpiexec,
  args: ['-np', '2', test_mumps_exe, files('input_simpletest_real.txt')],
  is_parallel: false,
  env: nomalloc_env,
  suite: 'unit')

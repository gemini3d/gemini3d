pde_elliptic = library('pde_elliptic', 'PDEelliptic.F90',
  dependencies: [mpi, mumps],
  link_with: const,
  fortran_args: real_bits)

potential = library('potential',
  sources: ['potential_comm_mumps.f90', 'potential_mumps.F90', 'boundary_conditions/potentialBCs_mumps.f90'],
  link_with: [pde_elliptic, calculus, collisions, const, grid, interp, mpimod, timeutils],
  dependencies: [mumps, scalapack, blacs, mpi, threads],
  fortran_args: real_bits)

test_potential2d = executable('test_potential2d', 'test_potential2D.f90',
  link_with: [pde_elliptic, const],
  dependencies: [blacs, mpi],
  link_language: 'fortran')
test('potential2d', test_potential2d, timeout: 30, suite: 'unit')

test('potential2d_octave', octave,
  args: ['-q', '--eval', 'test_potential2D("' + meson.build_root() / 'test_potential2D.dat")'],
  timeout: 15,
  workdir: meson.current_source_dir(),
  suite: 'unit', priority: -100)

test('potential2d_python', python,
  args: ['test_potential2d.py', meson.build_root() / 'test_potential2D.dat'],
  timeout: 15,
  workdir: meson.current_source_dir(),
  suite: 'unit', priority: -100)

#test_potential3d = executable('test_potential3d', 'test_potential3D.f90',
#  dependencies: [mpi, mumps])
#test('potential3d', test_potential3d, timeout: 60, suite: 'unit')

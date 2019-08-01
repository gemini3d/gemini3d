mpimod = library('mpimod',
  'mpimod/mpimod.F90', 'mpimod/mpisend.f90', 'mpimod/mpirecv.f90', 'mpimod/mpihalo.f90',
  dependencies: mpi,
  link_with: const,
  fortran_args: real_bits)

grid = library('grid', 'grid/grid.f90',
  dependencies: mpi,
  link_with: [const, mpimod])

calculus = library('calculus', 'calculus/calculus.f90', 'interpolation/interpolation.f90',
  link_with: [const, grid, mpimod],
  dependencies : mpi,
  fortran_args: real_bits)

# FIXME: debug Meson internals
# mumps = dependency('MUMPS', cmake_module_path : 'cmake/Modules')
mumps_root = get_option('MUMPS_ROOT')
scalapack_root = get_option('SCALAPACK_ROOT')

if mumps_root == ''
  mumpsinc = '/usr/include'
  mumpslib = [fc.find_library('dmumps', required : true)]
else
  message('searching for MUMPS in ' + mumps_root)
  mumpsinc = mumps_root / 'include'
  mumpslib = [fc.find_library('dmumps', required : true, dirs : mumps_root),
              fc.find_library('mumps_common', required : true, dirs : mumps_root)]
endif
mumps = declare_dependency(include_directories: mumpsinc, dependencies: mumpslib)

if scalapack_root == ''
  scalapack = fc.find_library('scalapack-openmpi', required: false)
  if not scalapack.found()
    scalapack = fc.find_library('scalapack')
  endif
else
  message('searching for SCALAPACK in ' + scalapack_root)
  scalapack = fc.find_library('scalapack', dirs : scalapack_root / 'lib')
endif

blacs = fc.find_library('blacs-openmpi', required: false)

threads = dependency('threads')

potential = library('potential', 'potential/potential_comm_mumps.f90', 'potential/potential_mumps.F90', 'potential/boundary_conditions/potentialBCs_mumps.f90',
  link_with: [calculus, collisions, fsutils, mpimod, const, grid, timeutils],
  dependencies: [mumps, scalapack, blacs, mpi, threads],
  fortran_args: real_bits)

project('GEMINI', 'fortran',
  default_options : ['default_library=static'])
  
REALBITS = '-DREALBITS='+get_option('realbits')

fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_global_arguments('-O3','-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function', language : 'fortran')
endif

mpi = dependency('mpi', language : 'fortran', required : true)
#mumps = dependency('mumps')#s, module_path : 'cmake/Modules')
mumps = fc.find_library('dmumps', required : true)
lapack = dependency('lapack', required : true)


calculus = library('calculus', 'numerical/constants/phys_consts.F90', 'numerical/calculus/calculus.f90', 'numerical/grid/grid.f90', 'numerical/mpimod/mpimod.F90', 'numerical/interpolation/interpolation.f90', 'vendor/lapack95/gbsv.F90',
  dependencies : mpi,
  fortran_args: REALBITS)

io = library('io','io/io.F90', 'io/expanduser.f90',
  link_with: calculus,
  fortran_args: REALBITS)
  
msis = library('msis', 'vendor/msis00/msis00_gfortran.f',
  fortran_args : '-w')

neutral = library('neutral', 'neutral/neutral.f90', 'temporal/timeutils.f90',
  dependencies : mpi,
  link_with : [calculus, io, msis])

gemini = executable('gemini', 'gemini.f90', 'temporal/temporal.f90', 'numerical/potential/potential_comm_mumps.f90', 'numerical/potential/potential_mumps.F90', 'collisions/collisions.f90', 'numerical/potential/boundary_conditions/potentialBCs_mumps.f90', 'multifluid/multifluid.f90', 
'ionization/ionization.f90', 'ionization/boundary_conditions/precipBCs_mod.f90', 'ionization/glow_dummy.f90',
'sources/sources.f90', 'numerical/advection/advec_mpi.f90', 'numerical/diffusion/diffusion.F90',

  link_with : [io, neutral, calculus],
  dependencies : [mpi, mumps, lapack],
  include_directories: '/usr/include',
  fortran_args: REALBITS)

test('gemini 2D', gemini,
  args : ['../initialize/test2d/config.ini', '/tmp/test2d'])

#hdf5 = dependency('hdf5')

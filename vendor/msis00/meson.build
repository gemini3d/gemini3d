if fc.get_id() == 'gcc'
  oldflags = ['-std=legacy', quiet, '-fno-implicit-none',
              '-Wno-error=array-bounds', '-fcheck=no-all']
elif fc.get_id() == 'intel'
  oldflags = [quiet, '-nocheck']
elif fc.get_id() == 'intel-cl'
  oldflags = [quiet, '/nocheck']
else
  oldflags = [quiet]
endif

msis = library('msis', 'msis00_gfortran.f',
  fortran_args : oldflags)

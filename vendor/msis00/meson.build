
if fcid == 'gcc'
  oldflags = ['-std=legacy', '-w', '-fno-implicit-none',
              '-Wno-error=array-bounds', '-fcheck=no-all']
elif fcid == 'intel'
  oldflags = ['-nowarn', '-nocheck']
elif fcid == 'intel-cl'
  oldflags = ['/nowarn', '/nocheck']
else
  oldflags = ['-w']
endif

msis = library('msis', 'msis00_gfortran.f',
  fortran_args : oldflags)

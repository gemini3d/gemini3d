
if fcid == 'gcc'
  oldflags = ['-std=legacy', '-w', '-fno-implicit-none']
elif fcid == 'intel'
  oldflags = '-nowarn'
elif fcid == 'intel-cl'
  oldflags = '/nowarn'
else
  oldflags = ['-w']
endif

msis = library('msis', 'msis00_gfortran.f',
  fortran_args : oldflags)

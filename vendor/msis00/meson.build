
if fc.get_id() == 'gcc'
  oldflags = ['-std=legacy', '-w', '-fno-implicit-none']
elif fc.get_id() == 'intel'
  oldflags = '-nowarn'
elif fc.get_id() == 'intel-cl'
  oldflags == '/nowarn'
else
  oldflags = ['-w']
endif

msis = library('msis', 'msis00_gfortran.f',
  fortran_args : oldflags)

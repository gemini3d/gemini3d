io_src = files(
'io.f90',
'aurora.f90',
'mag.f90',
'plasma.f90',
'output.f90',
'input.f90'
)

# accomodate intel compiler non-standard inquire(directory=)
if ['intel', 'intel-cl'].contains(fc.get_id())
  io_src += files('path_exists_intel.f90')
else
  io_src += files('path_exists.f90')
endif

# FIXME: when validated, use if hdf5_interface.found() as for output.
io_src += files('plasma_input_raw.f90')

if hdf5_interface.found()
  io_src += files('plasma_output_hdf5.f90', 'mag_hdf5.f90', 'aurora_hdf5.f90')
elif netcdf_interface.found()
  io_src += files('plasma_output_ncdf.f90', 'mag_ncdf.f90', 'aurora_ncdf.f90')
else
  io_src += files('plasma_output_raw.f90', 'mag_raw.f90', 'aurora_raw.f90')
endif

fsutils = library('fsutils',
  sources: ['mkdir.F90', 'expanduser.f90'])

logging = library('logging', 'logging.f90')

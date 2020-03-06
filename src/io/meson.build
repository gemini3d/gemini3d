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

if hdf5_interface.found()
  io_src += files('plasma_input_hdf5.f90', 'plasma_output_hdf5.f90', 'mag_hdf5.f90', 'aurora_hdf5.f90')
elif netcdf_interface.found()
  io_src += files('plasma_input_ncdf.f90', 'plasma_output_ncdf.f90', 'mag_ncdf.f90', 'aurora_ncdf.f90')
else
  io_src += files('plasma_input_raw.f90', 'plasma_output_raw.f90', 'mag_raw.f90', 'aurora_raw.f90')
endif

# === reader

if hdf5_interface.found()
  reader_src = files('reader_hdf5.f90')
  reader_dep = hdf5_interface
elif netcdf_interface.found()
  reader_src = files('reader_nc4.f90')
  reader_dep = netcdf_interface
else
  reader_src = files('reader_raw.f90')
  reader_dep = []
endif
reader = library('reader',
  sources: ['reader.f90', reader_src],
  link_with: [const],
  dependencies: reader_dep)

# == other

pathlib = library('pathlib', 'pathlib.F90')

logging = library('logging', 'logging.f90')

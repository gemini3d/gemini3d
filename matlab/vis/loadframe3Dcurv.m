function dat = loadframe3Dcurv(filename)

narginchk(1,1)
[~,~,ext] = fileparts(filename);

switch ext
  case '.h5', dat = loadframe3Dcurv_hdf5(filename);
  case '.dat', dat = loadframe3Dcurv_raw(filename);
  case '.nc', error('NetCDF4 not yet handled in Matlab')
  otherwise, error(['unknown file type', filename])
end

end % function

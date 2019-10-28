function [ne,v1,Ti,Te,J1,v2,v3,J2,J3,ns,vs1,Ts,Phitop] = loadframe3Dcurv(direc, filename)

[~,~,ext] = fileparts(filename);

switch ext
  case '.dat', [ne,v1,Ti,Te,J1,v2,v3,J2,J3,ns,vs1,Ts,Phitop] = loadframe3Dcurv_raw(direc, filename);
  case '.h5', [ne,v1,Ti,Te,J1,v2,v3,J2,J3,ns,vs1,Ts,Phitop] = loadframe3Dcurv_hdf5(direc, filename);
  otherwise, error(['unknown file type', filename])
end

end % function

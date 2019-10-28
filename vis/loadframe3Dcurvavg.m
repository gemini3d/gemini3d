function [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg(direc, filename)

narginchk(2,2)
[~,~,ext] = fileparts(filename);

switch ext
  case '.dat', [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg_raw(direc, filename);
  case '.h5', [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg_hdf5(direc, filename);
  otherwise, error(['unknown file type ',filename])
end

end % function

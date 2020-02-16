function [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg(direc, filename)

narginchk(1,2)
if nargin == 2
  filename = [direc, filesep, filename];
else
  filename = direc;
end

[~,~,ext] = fileparts(filename);

switch ext
  case '.dat', [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg_raw(filename);
  case '.h5', [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg_hdf5(filename);
  otherwise, error(['unknown file type ',filename])
end

end % function

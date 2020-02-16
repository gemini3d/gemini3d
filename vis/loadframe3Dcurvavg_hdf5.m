function [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg_hdf5(filename)

narginchk(1,1)
%% SIMULATION SIZE
lxs = simsize(filename);
%% SIMULATION GRID FILE
% (NOTE THAT THIS IS NOT THE ENTIRE THING - THAT NEEDS TO BE DONE WITH READGRID.M.  WE NEED THIS HERE TO DO MESHGRIDS
% [x1, x2, x3] = simaxes(filename);
%% SIMULATIONS RESULTS
assert(is_file(filename), [filename,' does not exist '])

% simdate=zeros(1,6);    %datevec-style array

if isoctave
  D = load(filename);
%   simdate(1:3) = D.time.ymd;
%   simdate(4) = D.time.UThour;
  ne = D.neall;
  v1 = D.v1avgall;
  Ti = D.Tavgall;
  Te = D.TEall;
  J1 = D.J1all;
  J2 = D.J2all;
  J3 = D.J3all;
  v2 = D.v2avgall;
  v3 = D.v3avgall;
  Phitop = D.Phiall;
else
%   simdate(1:3) = h5read(filename, '/time/ymd');
%   simdate(4) = h5read(filename, '/time/UThour');
  %% Number densities
  ne = h5read(filename, '/neall');
  %% Parallel Velocities
  v1 = h5read(filename, '/v1avgall');
  %% Temperatures
  Ti = h5read(filename, '/Tavgall');
  Te = h5read(filename, '/TEall');
  %% Current densities
  J1 = h5read(filename, '/J1all');
  J2 = permute(h5read(filename, '/J2all'), [1,3,2]);
  J3 = permute(h5read(filename, '/J3all'), [1,3,2]);
  %% Perpendicular drifts
  v2 = h5read(filename, '/v2avgall');
  v3 = h5read(filename, '/v3avgall');
  %% Topside potential
  Phitop = h5read(filename, '/Phiall');
end
%% REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if any(lxs(2:3) == 1)    %a 2D simulations was done in x1 and x3
  ne = squeeze(ne);
  v1 = squeeze(v1);
  Ti = squeeze(Ti);
  Te = squeeze(Te);
  J1 = squeeze(J1);
  J2 = squeeze(J2);
  J3 = squeeze(J3);
  v2 = squeeze(v2);
  v3 = squeeze(v3);
end

end % function

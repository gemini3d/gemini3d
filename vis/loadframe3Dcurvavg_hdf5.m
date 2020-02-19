function dat = loadframe3Dcurvavg_hdf5(filename)

narginchk(1,1)
%% SIMULATION SIZE
lxs = simsize(filename);
%% SIMULATION GRID FILE
% (NOTE THAT THIS IS NOT THE ENTIRE THING - THAT NEEDS TO BE DONE WITH READGRID.M.  WE NEED THIS HERE TO DO MESHGRIDS
% [x1, x2, x3] = simaxes(filename);
%% SIMULATIONS RESULTS
assert(is_file(filename), [filename,' does not exist '])
dat.filename = filename;

% simdate=zeros(1,6);    %datevec-style array

if isoctave
  D = load(filename);
%   simdate(1:3) = D.time.ymd;
%   simdate(4) = D.time.UThour;
  dat.ne = D.neall;
  dat.v1 = D.v1avgall;
  dat.Ti = D.Tavgall;
  dat.Te = D.TEall;
  dat.J1 = D.J1all;
  dat.J2 = D.J2all;
  dat.J3 = D.J3all;
  dat.v2 = D.v2avgall;
  dat.v3 = D.v3avgall;
  dat.Phitop = D.Phiall;
else
%   simdate(1:3) = h5read(filename, '/time/ymd');
%   simdate(4) = h5read(filename, '/time/UThour');
  %% Number densities
  dat.ne = h5read(filename, '/neall');
  %% Parallel Velocities
  dat.v1 = h5read(filename, '/v1avgall');
  %% Temperatures
  dat.Ti = h5read(filename, '/Tavgall');
  dat.Te = h5read(filename, '/TEall');
  %% Current densities
  dat.J1 = h5read(filename, '/J1all');
  dat.J2 = permute(h5read(filename, '/J2all'), [1,3,2]);
  dat.J3 = permute(h5read(filename, '/J3all'), [1,3,2]);
  %% Perpendicular drifts
  dat.v2 = h5read(filename, '/v2avgall');
  dat.v3 = h5read(filename, '/v3avgall');
  %% Topside potential
  dat.Phitop = h5read(filename, '/Phiall');
end
%% REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if any(lxs(2:3) == 1)    %a 2D simulations was done in x1 and x3
  dat.ne = squeeze(dat.ne);
  dat.v1 = squeeze(dat.v1);
  dat.Ti = squeeze(dat.Ti);
  dat.Te = squeeze(dat.Te);
  dat.J1 = squeeze(dat.J1);
  dat.J2 = squeeze(dat.J2);
  dat.J3 = squeeze(dat.J3);
  dat.v2 = squeeze(dat.v2);
  dat.v3 = squeeze(dat.v3);
end

end % function

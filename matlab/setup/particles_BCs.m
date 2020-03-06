function particles_BCs(p, xg)
% create particle precipitation
narginchk(2,2)
validateattributes(p, {'struct'}, {'scalar'})
validateattributes(xg, {'struct'}, {'scalar'})
%% CREATE PRECIPITATION CHARACTERISTICS data
% number of grid cells.
% This will be interpolated to grid, so 100x100 is arbitrary
llon=100;
llat=100;

if xg.lx(2) == 1    % cartesian
  llon=1;
elseif xg.lx(3) == 1
  llat=1;
end

%% TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
% dtprec is set in config.nml
time = 0:p.dtprec:p.tdur;
Nt = numel(time);

%% time
UTsec = p.UTsec0 + time;     % seconds from beginning of hour
UThrs = UTsec/3600;
expdate=cat(2,repmat(p.ymd(:)',[Nt,1]),UThrs(:),zeros(Nt,1),zeros(Nt,1));

%% CREATE PRECIPITATION INPUT DATA
% Qit: energy flux [mW m^-2]
% E0it: characteristic energy [eV]
Qit = zeros(llon, llat, Nt);
E0it = zeros(llon,llat, Nt);

% did user specify on/off time? if not, assume always on.
if isfield(p, 'precip_startsec')
  [~, i_on] = min(abs(time - p.precip_startsec));
else
  i_on = 1;
end

if isfield(p, 'precip_endsec')
  [~, i_off] = min(abs(time - p.precip_endsec));
else
  i_off = Nt;
end

pg = precip_grid(xg, p, llat, llon);

for i = i_on:i_off
   Qit(:,:,i) = precip_gaussian2d(pg);
  E0it(:,:,i) = 5e3;
end

%% CONVERT THE ENERGY TO EV
%E0it = max(E0it,0.100);
%E0it = E0it*1e3;

%% SAVE to files
% LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
% FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.
% THE EFIELD DATA DO NOT NEED TO BE SMOOTHED.

outdir = absolute_path([p.simdir, '/inputs/prec_inputs/']);
makedir(outdir)

switch p.format
  case {'h5','hdf5'}, write_hdf5(outdir, llon, llat, pg.mlon, pg.mlat, expdate, Nt, Qit, E0it)
  case {'dat','raw'}, write_raw(outdir, llon, llat, pg.mlon, pg.mlat, expdate, Nt, Qit, E0it, p.realbits)
  otherwise, error(['unknown file format ', p.format])
end

end % function


function write_hdf5(outdir, llon, llat, mlon, mlat, expdate, Nt, Qit, E0it)
narginchk(10,10)

fn = [outdir, '/simsize.h5'];
disp(['write ', fn])
h5save(fn, '/llon', int32(llon))
h5save(fn, '/llat', int32(llat))

freal = 'float32';

fn = [outdir, '/simgrid.h5'];
disp(['write ', fn])
h5save(fn, '/mlon', mlon, [], freal)
h5save(fn, '/mlat', mlat, [], freal)

for i = 1:Nt
  UTsec = expdate(i,4)*3600 + expdate(i,5)*60 + expdate(i,6);
  ymd = expdate(i, 1:3);

  fn = [outdir, filesep, datelab(ymd,UTsec), '.h5'];
  disp(['writing ', fn])

  h5save(fn, '/Qp', Qit(:,:,i), [], freal)
  h5save(fn, '/E0p', E0it(:,:,i), [], freal)
end

end % function


function write_raw(outdir, llon, llat, mlon, mlat, expdate, Nt, Qit, E0it, realbits)
narginchk(10,10)

filename=[outdir, '/simsize.dat'];
disp(['write ', filename])
fid=fopen(filename, 'w');
fwrite(fid,llon,'integer*4');
fwrite(fid,llat,'integer*4');
fclose(fid);

freal = ['float', int2str(realbits)];

filename=[outdir, '/simgrid.dat'];
disp(['write ', filename])

fid=fopen(filename,'w');
fwrite(fid,mlon, freal);
fwrite(fid,mlat, freal);
fclose(fid);

for i = 1:Nt
  UTsec = expdate(i,4)*3600 + expdate(i,5)*60 + expdate(i,6);
  ymd = expdate(i, 1:3);

  filename = [outdir, filesep, datelab(ymd,UTsec), '.dat'];
  disp(['writing ', filename])

  fid = fopen(filename,'w');
  fwrite(fid,Qit(:,:,i), freal);
  fwrite(fid,E0it(:,:,i), freal);
  fclose(fid);
end

end % function

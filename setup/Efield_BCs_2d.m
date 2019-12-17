function E = Efield_BCs_2d(dir_grid, format, dir_config, realbits)

narginchk(2, 4)
validateattributes(dir_grid, {'char'}, {'vector'})
validateattributes(format, {'char'}, {'vector'})

cwd = fileparts(mfilename('fullpath'));
if nargin < 3 || isempty(dir_config), dir_config = cwd; end
validateattributes(dir_config, {'char'}, {'vector'})

if nargin < 4, realbits = 64; end
validateattributes(realbits, {'numeric'}, {'scalar', 'integer', 'positive'})

dir_grid = absolute_path(dir_grid);
dir_out = [dir_grid, '/Efield_inputs'];

if ~isfolder(dir_out)
  mkdir(dir_out);
end


%% READ IN THE SIMULATION INFORMATION
[ymd0, UTsec0] = readconfig(dir_config);

xg = readgrid(dir_grid);
lx1 = xg.lx(1);
lx2 = xg.lx(2);
lx3 = xg.lx(3);


%CREATE A 'DATASET' OF ELECTRIC FIELD INFO
llon=100;
llat=100;
if xg.lx(2) == 1    %this is cartesian-specific code
  llon=1;
elseif xg.lx(3) == 1
  llat=1;
end
thetamin=min(xg.theta(:));
thetamax=max(xg.theta(:));
mlatmin=90-thetamax*180/pi;
mlatmax=90-thetamin*180/pi;
mlonmin=min(xg.phi(:))*180/pi;
mlonmax=max(xg.phi(:))*180/pi;
latbuf=1/100*(mlatmax-mlatmin);
lonbuf=1/100*(mlonmax-mlonmin);
E.mlat = linspace(mlatmin-latbuf,mlatmax+latbuf,llat);
E.mlon = linspace(mlonmin-lonbuf,mlonmax+lonbuf,llon);
[E.MLON, E.MLAT] = ndgrid(E.mlon, E.mlat);
mlonmean = mean(E.mlon);
% mlatmean=mean(E.mlat);

%% WIDTH OF THE DISTURBANCE
fracwidth = 1/7;
% mlatsig = fracwidth*(mlatmax-mlatmin);
mlonsig=fracwidth*(mlonmax-mlonmin);
sigx2=fracwidth*(max(xg.x2)-min(xg.x2));

%% TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
tmin=0;
tmax=300;
dt = 1.;  % [seconds]  % FIXME: shouldn't this be from config.nml
time = tmin:dt:tmax;
lt = length(time);
%% SET UP TIME VARIABLES
ymd = ymd0;
UTsec = UTsec0+time;     %time given in file is the seconds from beginning of hour
UThrs = UTsec/3600;
E.expdate = cat(2, repmat(ymd(:)',[lt,1]), UThrs', zeros(lt,1), zeros(lt,1));
% t=datenum(E.expdate);

%% CREATE DATA FOR BACKGROUND ELECTRIC FIELDS
E.Exit=zeros(llon,llat,lt);
E.Eyit=zeros(llon,llat,lt);
for it=1:lt
  E.Exit(:,:,it)=zeros(llon,llat);   %V/m
  E.Eyit(:,:,it)=zeros(llon,llat);
end

%% CREATE DATA FOR BOUNDARY CONDITIONS FOR POTENTIAL SOLUTION
flagdirich=1;   %if 0 data is interpreted as FAC, else we interpret it as potential
E.Vminx1it=zeros(llon,llat,lt);
E.Vmaxx1it=zeros(llon,llat,lt);
E.Vminx2ist=zeros(llat,lt);
E.Vmaxx2ist=zeros(llat,lt);
E.Vminx3ist=zeros(llon,lt);
E.Vmaxx3ist=zeros(llon,lt);
Etarg=50e-3;            % target E value in V/m
if lx3 == 1
  pk = Etarg*sigx2 .* xg.h2(lx1, floor(lx2/2), 1).*sqrt(pi)./2;
elseif lx2 == 1
  pk = Etarg*sigx2 .* xg.h2(lx1, 1, floor(lx3/2)).*sqrt(pi)./2;
end
% x2ctr = 1/2*(xg.x2(lx2)+xg.x2(1));
for it=1:lt
  E.Vminx1it(:,:,it)=zeros(llon,llat);
  E.Vmaxx1it(:,:,it)=pk.*erf((E.MLON-mlonmean)/mlonsig);%.*erf((MLAT-mlatmean)/mlatsig);
  E.Vminx2ist(:,it)=zeros(llat,1);     %these are just slices
  E.Vmaxx2ist(:,it)=zeros(llat,1);
  E.Vminx3ist(:,it)=zeros(llon,1);
  E.Vmaxx3ist(:,it)=zeros(llon,1);
end

%% check for NaNs
% this is also done in Fortran, but just to help ensure results.
assert(all(isfinite(E.Exit(:))), 'NaN in Exit')
assert(all(isfinite(E.Eyit(:))), 'NaN in Eyit')
assert(all(isfinite(E.Vminx1it(:))), 'NaN in Vminxlit')
assert(all(isfinite(E.Vmaxx1it(:))), 'NaN in Vmaxxlit')
assert(all(isfinite(E.Vminx2it(:))), 'NaN in Vminx2it')
assert(all(isfinite(E.Vmaxx2it(:))), 'NaN in Vmaxx2it')
assert(all(isfinite(E.Vminx3it(:))), 'NaN in Vminx3it')
assert(all(isfinite(E.Vmaxx3it(:))), 'NaN in Vmaxx3it')
%% SAVE THESE DATA TO APPROPRIATE FILES
% LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
% FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.
% THE EFIELD DATA DO NOT TYPICALLY NEED TO BE SMOOTHED.

switch format
  case {'raw', 'dat'}, writeraw(dir_out, realbits)
  case {'h5', 'hdf5'}, writehdf5(dir_out)
  otherwise, error(['unknown data format ', format])
end


if ~nargout, clear('E'), end

%% nested functions

function writehdf5(dir_out)

fn = [dir_out, '/simsize.h5'];
if isfile(fn), delete(fn), end
h5save(fn, '/llon', llon)
h5save(fn, '/llat', llat)

fn = [dir_out, '/simgrid.h5'];
if isfile(fn), delete(fn), end
h5save(fn, '/mlon', E.mlon)
h5save(fn, '/mlat', E.mlat)

for i = 1:lt
  UTsec = E.expdate(i, 4)*3600 + E.expdate(i,5)*60 + E.expdate(i,6);
  ymd = E.expdate(i, 1:3);

  fn = [dir_out, filesep, datelab(ymd,UTsec), '.h5'];
  disp(['write: ', fn])

  %FOR EACH FRAME WRITE A BC TYPE AND THEN OUTPUT BACKGROUND AND BCs
  h5save(fn, '/flagdirich', flagdirich)
  h5save(fn, '/Exit', E.Exit(:,:,i))
  h5save(fn, '/Eyit', E.Eyit(:,:,i))
  h5save(fn, '/Vminx1it', E.Vminx1it(:,:,i))
  h5save(fn, '/Vmaxx1it', E.Vmaxx1it(:,:,i))
  h5save(fn, '/Vminx2ist', E.Vminx2ist(:,i))
  h5save(fn, '/Vmaxx2ist', E.Vmaxx2ist(:,i))
  h5save(fn, '/Vminx3ist', E.Vminx3ist(:,i))
  h5save(fn, '/Vmaxx3ist', E.Vmaxx3ist(:,i))
end
end % function


function writeraw(dir_out, realbits)

assert(any(realbits == [32, 64]), 'realbits == 32 or 64')

freal = ['float',int2str(realbits)];

fid = fopen([dir_out, '/simsize.dat'], 'w');
fwrite(fid, llon, 'integer*4');
fwrite(fid, llat, 'integer*4');
fclose(fid);

fid=fopen([dir_out, '/simgrid.dat'], 'w');
fwrite(fid, E.mlon, freal);
fwrite(fid, E.mlat, freal);
fclose(fid);

for i=1:lt
  UTsec = E.expdate(i,4)*3600 + E.expdate(i,5)*60 + E.expdate(i,6);
  ymd = E.expdate(i,1:3);
  filename = [dir_out, filesep, datelab(ymd,UTsec), '.dat'];
  disp(['write: ',filename])
  fid = fopen(filename, 'w');

  %FOR EACH FRAME WRITE A BC TYPE AND THEN OUTPUT BACKGROUND AND BCs
  fwrite(fid, flagdirich, freal);
  fwrite(fid, E.Exit(:,:,i), freal);
  fwrite(fid, E.Eyit(:,:,i), freal);
  fwrite(fid, E.Vminx1it(:,:,i), freal);
  fwrite(fid, E.Vmaxx1it(:,:,i), freal);
  fwrite(fid, E.Vminx2ist(:,i), freal);
  fwrite(fid, E.Vmaxx2ist(:,i), freal);
  fwrite(fid, E.Vminx3ist(:,i), freal);
  fwrite(fid, E.Vmaxx3ist(:,i), freal);

  fclose(fid);
end

end % function


end % function

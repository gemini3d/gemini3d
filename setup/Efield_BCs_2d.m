function E = Efield_BCs_2d(p, xg)

narginchk(2, 2)
validateattributes(p, {'struct'}, {'scalar'}, mfilename, 'sim parameters', 1)
validateattributes(xg, {'struct'}, {'scalar'})

dir_out = absolute_path([p.simdir, '/inputs/Efield_inputs']);
makedir(dir_out);

lx1 = xg.lx(1);
lx2 = xg.lx(2);
lx3 = xg.lx(3);

%% CREATE ELECTRIC FIELD DATASET
E.llon=100;
E.llat=100;
% NOTE: cartesian-specific code
if lx2 == 1
  E.llon = 1;
elseif lx3 == 1
  E.llat = 1;
end
thetamin = min(xg.theta(:));
thetamax = max(xg.theta(:));
mlatmin = 90-thetamax*180/pi;
mlatmax = 90-thetamin*180/pi;
mlonmin = min(xg.phi(:))*180/pi;
mlonmax = max(xg.phi(:))*180/pi;

% add a 1% buff
latbuf = 1/100 * (mlatmax-mlatmin);
lonbuf = 1/100 * (mlonmax-mlonmin);
E.mlat = linspace(mlatmin-latbuf, mlatmax+latbuf, E.llat);
E.mlon = linspace(mlonmin-lonbuf, mlonmax+lonbuf, E.llon);
[E.MLON, E.MLAT] = ndgrid(E.mlon, E.mlat);
mlonmean = mean(E.mlon);
mlatmean = mean(E.mlat);

%% WIDTH OF THE DISTURBANCE
mlatsig = p.Efield_fracwidth*(mlatmax-mlatmin);
mlonsig = p.Efield_fracwidth*(mlonmax-mlonmin);
sigx2 = p.Efield_fracwidth*(max(xg.x2)-min(xg.x2));
sigx3 = p.Efield_fracwidth*(max(xg.x3)-min(xg.x3));
%% TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
tmin = 0;
time = tmin:p.dtE0:p.tdur;
Nt = length(time);
%% SET UP TIME VARIABLES
UTsec = p.UTsec0 + time;     %time given in file is the seconds from beginning of hour
UThrs = UTsec / 3600;
E.expdate = cat(2, repmat(p.ymd(:)',[Nt, 1]), UThrs', zeros(Nt, 1), zeros(Nt, 1));
% t=datenum(E.expdate);
%% CREATE DATA FOR BACKGROUND ELECTRIC FIELDS
E.Exit = zeros(E.llon, E.llat, Nt);
E.Eyit = zeros(E.llon, E.llat, Nt);

% put custom E-field background in this for loop
%{
for it=1:Nt
  E.Exit(:,:,it) = ;   %V/m
  E.Eyit(:,:,it) = ;
end
%}
%% CREATE DATA FOR BOUNDARY CONDITIONS FOR POTENTIAL SOLUTION
p.flagdirich = 1;   %if 0 data is interpreted as FAC, else we interpret it as potential
E.Vminx1it = zeros(E.llon,E.llat, Nt);
E.Vmaxx1it = zeros(E.llon,E.llat, Nt);
%these are just slices
E.Vminx2ist = zeros(E.llat, Nt);
E.Vmaxx2ist = zeros(E.llat, Nt);
E.Vminx3ist = zeros(E.llon, Nt);
E.Vmaxx3ist = zeros(E.llon, Nt);

if p.Etarg > 1
  warning(['Etarg units V/m but is ',num2str(p.Etarg), ' V/m realistic?'])
end

if lx3 == 1 % east-west
  pk = p.Etarg * sigx2 .* xg.h2(lx1, floor(lx2/2), 1) .* sqrt(pi)./2;
elseif lx2 == 1 % north-south
  pk = p.Etarg * sigx3 .* xg.h3(lx1, 1, floor(lx3/2)) .* sqrt(pi)./2;
end

% x2ctr = 1/2*(xg.x2(lx2) + xg.x2(1));
for i = 1:Nt
  % put your functions in these if you want
  %{
  E.Vminx1it(:,:,i) = ;
  %}
  if lx2 == 1
    E.Vmaxx1it(:,:,i) = pk .* erf((E.MLAT - mlatmean)/mlatsig);
  elseif lx3 == 1
    E.Vmaxx1it(:,:,i) = pk .* erf((E.MLON - mlonmean)/mlonsig);
  end
  % put your functions in these if you want
  %{
  E.Vminx2ist(:,i) = ;
  E.Vmaxx2ist(:,i) = ;
  E.Vminx3ist(:,i) = ;
  E.Vmaxx3ist(:,i) = ;
  %}
end

%% check for NaNs
% this is also done in Fortran, but just to help ensure results.
assert(all(isfinite(E.Exit(:))), 'NaN in Exit')
assert(all(isfinite(E.Eyit(:))), 'NaN in Eyit')
assert(all(isfinite(E.Vminx1it(:))), 'NaN in Vminx1it')
assert(all(isfinite(E.Vmaxx1it(:))), 'NaN in Vmaxx1it')
assert(all(isfinite(E.Vminx2ist(:))), 'NaN in Vminx2ist')
assert(all(isfinite(E.Vmaxx2ist(:))), 'NaN in Vmaxx2ist')
assert(all(isfinite(E.Vminx3ist(:))), 'NaN in Vminx3ist')
assert(all(isfinite(E.Vmaxx3ist(:))), 'NaN in Vmaxx3ist')

%% SAVE THESE DATA TO APPROPRIATE FILES
% LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
% FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.
% THE EFIELD DATA DO NOT TYPICALLY NEED TO BE SMOOTHED.

switch p.format
  case {'raw', 'dat'}, writeraw(dir_out, E, p)
  case {'h5', 'hdf5'}, writehdf5(dir_out, E, p)
  otherwise, error(['unknown data format ', p.format])
end

if ~nargout, clear('E'), end

end % function


function writehdf5(dir_out, E, p)
narginchk(3,3)

fn = [dir_out, '/simsize.h5'];
if is_file(fn), delete(fn), end
h5save(fn, '/llon', int32(E.llon))
h5save(fn, '/llat', int32(E.llat))

fn = [dir_out, '/simgrid.h5'];
if is_file(fn), delete(fn), end

freal = 'float32';

h5save(fn, '/mlon', E.mlon, [], freal)
h5save(fn, '/mlat', E.mlat, [], freal)

disp(['write Efield to ', dir_out])
Nt = size(E.expdate, 1);
for i = 1:Nt
  UTsec = E.expdate(i, 4)*3600 + E.expdate(i,5)*60 + E.expdate(i,6);
  ymd = E.expdate(i, 1:3);

  fn = [dir_out, filesep, datelab(ymd,UTsec), '.h5'];

  %FOR EACH FRAME WRITE A BC TYPE AND THEN OUTPUT BACKGROUND AND BCs
  h5save(fn, '/flagdirich', int32(p.flagdirich))
  h5save(fn, '/Exit', E.Exit(:,:,i), [], freal)
  h5save(fn, '/Eyit', E.Eyit(:,:,i), [], freal)
  h5save(fn, '/Vminx1it', E.Vminx1it(:,:,i), [], freal)
  h5save(fn, '/Vmaxx1it', E.Vmaxx1it(:,:,i), [], freal)
  h5save(fn, '/Vminx2ist', E.Vminx2ist(:,i), [], freal)
  h5save(fn, '/Vmaxx2ist', E.Vmaxx2ist(:,i), [], freal)
  h5save(fn, '/Vminx3ist', E.Vminx3ist(:,i), [], freal)
  h5save(fn, '/Vmaxx3ist', E.Vmaxx3ist(:,i), [], freal)
end
end % function


function writeraw(dir_out, E, p)
narginchk(3,3)

freal = ['float', int2str(p.realbits)];

fid = fopen([dir_out, '/simsize.dat'], 'w');
fwrite(fid, E.llon, 'integer*4');
fwrite(fid, E.llat, 'integer*4');
fclose(fid);

fid = fopen([dir_out, '/simgrid.dat'], 'w');
fwrite(fid, E.mlon, freal);
fwrite(fid, E.mlat, freal);
fclose(fid);

Nt = size(E.expdate, 1);
for i = 1:Nt
  UTsec = E.expdate(i,4)*3600 + E.expdate(i,5)*60 + E.expdate(i,6);
  ymd = E.expdate(i,1:3);
  filename = [dir_out, filesep, datelab(ymd,UTsec), '.dat'];
  disp(['write: ',filename])
  fid = fopen(filename, 'w');

  %FOR EACH FRAME WRITE A BC TYPE AND THEN OUTPUT BACKGROUND AND BCs
%  fwrite(fid, p.flagdirich, 'int32');
  fwrite(fid, p.flagdirich, freal);   %FIXME - fortran code still wants this to be real...
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

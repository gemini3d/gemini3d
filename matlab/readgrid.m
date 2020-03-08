function xgf = readgrid(path, file_format, realbits)
%% READS A GRID FROM MATLAB
% OR POSSIBLY FORTRAN (THOUGH THIS IS NOT YET IMPLEMENTED AS OF 9/15/2016)
narginchk(1,3)
if nargin < 2 || isempty(file_format), file_format = 'auto'; end
validateattributes(file_format, {'char'}, {'vector'}, mfilename, 'raw or hdf5', 2)
if nargin < 3 || isempty(realbits), realbits = 64; end
validateattributes(realbits, {'numeric'}, {'scalar', 'integer'}, mfilename, '32 or 64', 3)

path = absolute_path(path);
if is_file(path)
  path = fileparts(path);
end

switch file_format
  case {'dat','raw'}, xgf = read_raw(path, realbits);
  case {'h5','hdf5'}, xgf = read_hdf5(path);
  otherwise
    if is_file([path, '/simsize.h5'])
      xgf = read_hdf5(path);
    elseif is_file([path,'/simsize.dat'])
      xgf = read_raw(path, realbits);
    else
      error(['no simsize file found in ',path])
    end
end

end % function


function xgf = read_hdf5(path)
for f = {[path, '/inputs/simsize.h5'], [path, '/simsize.h5']}
  sizefn = f{:};
  if is_file(sizefn)
    break
  end
end
assert(is_file(sizefn), [sizefn, ' not found'])

for f = {[path, '/inputs/simgrid.h5'], [path, '/simgrid.h5']}
  fn = f{:};
  if is_file(fn)
    break
  end
end
assert(is_file(fn), [fn, ' not found'])

xgf.filename = fn;

if isoctave
  L = load(sizefn);
  xgf = load(fn);

  try
    xgf.lx = L.lx;
  catch
    % octave bug: error: octave_base_value::int32_scalar_value(): wrong type argument 'int32 matrix'
    xgf.lx = [L.lx1; L.lx2; L.lx3];
  end

else

  try
    xgf.lx = h5read(sizefn, '/lx');
  catch
    xgf.lx = [h5read(sizefn, '/lx1'), h5read(sizefn, '/lx2'), h5read(sizefn, '/lx3')];
  end
  xgf.x1 = h5read(fn, '/x1');
  xgf.x1i = h5read(fn, '/x1i');
  xgf.dx1b = h5read(fn, '/dx1b');
  xgf.dx1h = h5read(fn, '/dx1h');
  xgf.x2 = h5read(fn, '/x2');
  xgf.x3 = h5read(fn, '/x3');

  xgf.h1 = h5read(fn, '/h1');
  xgf.h2 = h5read(fn, '/h2');
  xgf.h3 = h5read(fn, '/h3');

  xgf.alt = h5read(fn, '/alt');
  xgf.glat = h5read(fn, '/glat');
  xgf.glon = h5read(fn, '/glon');

  xgf.r = h5read(fn, '/r');
  xgf.theta = h5read(fn, '/theta');
  xgf.phi = h5read(fn, '/phi');
end

end  % function read_hdf5


function xgf = read_raw(path, realbits)


%% Size file
for f = {[path, '/inputs/simsize.dat'], [path, '/simsize.dat']}
  sizefn = f{:};
  if is_file(sizefn)
    break
  end
end
assert(is_file(sizefn), [sizefn, ' not found'])

%filename=[path, '/inputs/simsize.dat'];
filename=sizefn;
assert(is_file(filename), [filename,' is not a file.'])

fid=fopen(filename,'r');

xgf.filename = filename;
xgf.lx=fread(fid,3,'integer*4');
fclose(fid);
lx1=xgf.lx(1); lx2=xgf.lx(2); lx3=xgf.lx(3);
lgrid=lx1*lx2*lx3;
lgridghost=(lx1+4)*(lx2+4)*(lx3+4);
gridsize=[lx1,lx2,lx3];
gridsizeghost=[lx1+4,lx2+4,lx3+4];


%% Grid file
for f = {[path, '/inputs/simgrid.dat'], [path, '/simgrid.dat']}
  fn = f{:};
  if is_file(fn)
    break
  end
end
assert(is_file(fn), [fn, ' not found'])

%fin = [path, '/inputs/simgrid.dat'];
fin=fn;
assert(is_file(fin), [fin, ' is not a file.'])

freal = ['float',int2str(realbits)];

fid=fopen(fin,'r');

xgf.x1=fread(fid,lx1+4, freal);    %coordinate values
xgf.x1i=fread(fid,lx1+1, freal);
xgf.dx1b=fread(fid,lx1+3, freal);    %ZZZ - need to check that differences have appropriate ghost cell values, etc.
xgf.dx1h=fread(fid,lx1, freal);

xgf.x2=fread(fid,lx2+4, freal);
xgf.x2i=fread(fid,lx2+1, freal);
xgf.dx2b=fread(fid,lx2+3, freal);
xgf.dx2h=fread(fid,lx2, freal);

xgf.x3=fread(fid,lx3+4, freal);
xgf.x3i=fread(fid,lx3+1, freal);
xgf.dx3b=fread(fid,lx3+3, freal);
xgf.dx3h=fread(fid,lx3, freal);

tmp=fread(fid,lgridghost, freal);   %cell-centered metric coefficients
xgf.h1=reshape(tmp,gridsizeghost);
tmp=fread(fid,lgridghost, freal);
xgf.h2=reshape(tmp,gridsizeghost);
tmp=fread(fid,lgridghost, freal);
xgf.h3=reshape(tmp,gridsizeghost);

tmpsize=[lx1+1,lx2,lx3];
ltmp=prod(tmpsize);
tmp=fread(fid,ltmp, freal);    %interface metric coefficients
xgf.h1x1i=reshape(tmp,tmpsize);
tmp=fread(fid,ltmp, freal);
xgf.h2x1i=reshape(tmp,tmpsize);
tmp=fread(fid,ltmp, freal);
xgf.h3x1i=reshape(tmp,tmpsize);

tmpsize=[lx1,lx2+1,lx3];
ltmp=prod(tmpsize);
tmp=fread(fid,ltmp, freal);
xgf.h1x2i=reshape(tmp,tmpsize);
tmp=fread(fid,ltmp, freal);
xgf.h2x2i=reshape(tmp,tmpsize);
tmp=fread(fid,ltmp, freal);
xgf.h3x2i=reshape(tmp,tmpsize);

tmpsize=[lx1,lx2,lx3+1];
ltmp=prod(tmpsize);
tmp=fread(fid,ltmp, freal);
xgf.h1x3i=reshape(tmp,tmpsize);
tmp=fread(fid,ltmp, freal);
xgf.h2x3i=reshape(tmp,tmpsize);
tmp=fread(fid,ltmp, freal);
xgf.h3x3i=reshape(tmp,tmpsize);

%gravity, geographic coordinates, magnetic field strength? unit vectors?
tmp=fread(fid,lgrid, freal);    %gravitational field components
xgf.gx1=reshape(tmp,gridsize);
tmp=fread(fid,lgrid, freal);
xgf.gx2=reshape(tmp,gridsize);
tmp=fread(fid,lgrid, freal);
xgf.gx3=reshape(tmp,gridsize);

tmp=fread(fid,lgrid, freal);    %geographic coordinates
xgf.alt=reshape(tmp,gridsize);
tmp=fread(fid,lgrid, freal);
xgf.glat=reshape(tmp,gridsize);
tmp=fread(fid,lgrid, freal);
xgf.glon=reshape(tmp,gridsize);

tmp=fread(fid,lgrid, freal);    %magnetic field strength
xgf.Bmag=reshape(tmp,gridsize);

tmp=fread(fid,lx2*lx3, freal);    %magnetic field inclination - only one value for each field line
xgf.I=reshape(tmp,[lx2,lx3]);

tmp=fread(fid,lgrid, freal);    %points not to be solved
xgf.nullpts=reshape(tmp,gridsize);


%STUFF PAST THIS POINT ISN'T USED IN FORTRAN CODE BUT INCLUDED IN THE
%GRID FILE FOR COMPLETENESS
if ~feof(fid)
  tmpsize=[lx1,lx2,lx3,3];
  ltmp=prod(tmpsize);
  tmp=fread(fid,ltmp, freal);   %4D unit vectors (in cartesian components)

  if (feof(fid))    %for whatever reason, we sometimes don't hit eof until after first unit vector read...
     return;   %lazy as hell
  end

  xgf.e1=reshape(tmp,tmpsize);
  tmp=fread(fid,ltmp, freal);
  xgf.e2=reshape(tmp,tmpsize);
  tmp=fread(fid,ltmp, freal);
  xgf.e3=reshape(tmp,tmpsize);

  tmpsize=[lx1,lx2,lx3,3];
  ltmp=prod(tmpsize);
  tmp=fread(fid,ltmp, freal);
  xgf.er=reshape(tmp,tmpsize);
  tmp=fread(fid,ltmp, freal);
  xgf.etheta=reshape(tmp,tmpsize);
  tmp=fread(fid,ltmp, freal);
  xgf.ephi=reshape(tmp,tmpsize);

  tmp=fread(fid,lgrid, freal);    %spherical coordinates
  xgf.r=reshape(tmp,gridsize);
  tmp=fread(fid,lgrid, freal);
  xgf.theta=reshape(tmp,gridsize);
  tmp=fread(fid,lgrid, freal);
  xgf.phi=reshape(tmp,gridsize);

  tmp=fread(fid,lgrid, freal);     %cartesian coordinates
  xgf.x=reshape(tmp,gridsize);
  tmp=fread(fid,lgrid, freal);
  xgf.y=reshape(tmp,gridsize);
  tmp=fread(fid,lgrid, freal);
  xgf.z=reshape(tmp,gridsize);
end

fclose(fid);
end  % function read_raw

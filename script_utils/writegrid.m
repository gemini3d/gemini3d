function writegrid(p, xg)
%% write grid to raw binary files
% includes STUFF NOT NEEDED BY FORTRAN CODE BUT POSSIBLY USEFUL FOR PLOTTING

narginchk(2, 2)
validateattributes(p, {'struct'}, {'scalar'}, mfilename, 'simulation parameters', 1)
validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid parameters', 2)

outdir = absolute_path(p.simdir);
makedir(outdir)

switch p.format
  case {'dat','raw'}, write_raw(outdir, xg, p.realbits)
  case {'h5','hdf5'}, write_hdf5(outdir, xg, p.realbits)
  otherwise, error(['unknown file format ', p.format])
end

end % function


function write_hdf5(dir_out, xg, realbits)

fn = [dir_out, '/simsize.h5'];
disp(['write ',fn])
if is_file(fn), delete(fn), end
h5save(fn, '/lx1', int32(xg.lx(1)))
h5save(fn, '/lx2', int32(xg.lx(2)))
h5save(fn, '/lx3', int32(xg.lx(3)))

fn = [dir_out, '/simgrid.h5'];
disp(['write ',fn])
if is_file(fn), delete(fn), end

freal = ['float',int2str(realbits)];

h5save(fn, '/x1', xg.x1, [], freal)
h5save(fn, '/x1i', xg.x1i, [], freal)
h5save(fn, '/dx1b', xg.dx1b, [], freal)
h5save(fn, '/dx1h', xg.dx1h, [], freal)

h5save(fn, '/x2', xg.x2, [], freal)
h5save(fn, '/x2i', xg.x2i, [], freal)
h5save(fn, '/dx2b', xg.dx2b, [], freal)
h5save(fn, '/dx2h', xg.dx2h, [], freal)

h5save(fn, '/x3', xg.x3, [], freal)
h5save(fn, '/x3i', xg.x3i, [], freal)
h5save(fn, '/dx3b', xg.dx3b, [], freal)
h5save(fn, '/dx3h', xg.dx3h, [], freal)

h5save(fn, '/h1', xg.h1, [], freal)
h5save(fn, '/h2', xg.h2, [], freal)
h5save(fn, '/h3', xg.h3, [], freal)

h5save(fn, '/h1x1i', xg.h1x1i, [], freal)
h5save(fn, '/h2x1i', xg.h2x1i, [], freal)
h5save(fn, '/h3x1i', xg.h3x1i, [], freal)

h5save(fn, '/h1x2i', xg.h1x2i, [], freal)
h5save(fn, '/h2x2i', xg.h2x2i, [], freal)
h5save(fn, '/h3x2i', xg.h3x2i, [], freal)

h5save(fn, '/h1x3i', xg.h1x3i, [], freal)
h5save(fn, '/h2x3i', xg.h2x3i, [], freal)
h5save(fn, '/h3x3i', xg.h3x3i, [], freal)

h5save(fn, '/gx1', xg.gx1, [], freal)
h5save(fn, '/gx2', xg.gx2, [], freal)
h5save(fn, '/gx3', xg.gx3, [], freal)

h5save(fn, '/alt', xg.alt, [], freal)
h5save(fn, '/glat', xg.glat, [], freal)
h5save(fn, '/glon', xg.glon, [], freal)

h5save(fn, '/Bmag', xg.Bmag, [], freal)
h5save(fn, '/I', xg.I, [], freal)
h5save(fn, '/nullpts', xg.nullpts, [], freal)

h5save(fn, '/e1', xg.e1, [], freal)
h5save(fn, '/e2', xg.e2, [], freal)
h5save(fn, '/e3', xg.e3, [], freal)

h5save(fn, '/er', xg.er, [], freal)
h5save(fn, '/etheta', xg.etheta, [], freal)
h5save(fn, '/ephi', xg.ephi, [], freal)

h5save(fn, '/r', xg.r, [], freal)
h5save(fn, '/theta', xg.theta, [], freal)
h5save(fn, '/phi', xg.phi, [], freal)

h5save(fn, '/x', xg.x, [], freal)
h5save(fn, '/y', xg.y, [], freal)
h5save(fn, '/z', xg.z, [], freal)

end % function


function write_raw(outdir, xg, realbits)

freal = ['float', int2str(realbits)];

filename = [outdir, '/simsize.dat'];
disp(['write ',filename])
fid = fopen(filename, 'w');
fwrite(fid, xg.lx, 'integer*4');
fclose(fid);

fid = fopen([outdir, '/simgrid.dat'], 'w');

fwrite(fid,xg.x1, freal);    %coordinate values
fwrite(fid,xg.x1i, freal);
fwrite(fid,xg.dx1b, freal);
fwrite(fid,xg.dx1h, freal);

fwrite(fid,xg.x2, freal);
fwrite(fid,xg.x2i, freal);
fwrite(fid,xg.dx2b, freal);
fwrite(fid,xg.dx2h, freal);

fwrite(fid,xg.x3, freal);
fwrite(fid,xg.x3i, freal);
fwrite(fid,xg.dx3b, freal);
fwrite(fid,xg.dx3h, freal);

fwrite(fid,xg.h1, freal);   %cell-centered metric coefficients
fwrite(fid,xg.h2, freal);
fwrite(fid,xg.h3, freal);

fwrite(fid,xg.h1x1i, freal);    %interface metric coefficients
fwrite(fid,xg.h2x1i, freal);
fwrite(fid,xg.h3x1i, freal);

fwrite(fid,xg.h1x2i, freal);
fwrite(fid,xg.h2x2i, freal);
fwrite(fid,xg.h3x2i, freal);

fwrite(fid,xg.h1x3i, freal);
fwrite(fid,xg.h2x3i, freal);
fwrite(fid,xg.h3x3i, freal);

%gravity, geographic coordinates, magnetic field strength? unit vectors?
fwrite(fid,xg.gx1, freal);    %gravitational field components
fwrite(fid,xg.gx2, freal);
fwrite(fid,xg.gx3, freal);

fwrite(fid,xg.alt, freal);    %geographic coordinates
fwrite(fid,xg.glat, freal);
fwrite(fid,xg.glon, freal);

fwrite(fid,xg.Bmag, freal);    %magnetic field strength

fwrite(fid,xg.I, freal);    %magnetic field inclination

fwrite(fid,xg.nullpts, freal);    %points not to be solved


%NOT ALL OF THE REMAIN INFO IS USED IN THE FORTRAN CODE, BUT IT INCLUDED FOR COMPLETENESS
fwrite(fid,xg.e1, freal);   %4D unit vectors (in cartesian components)
fwrite(fid,xg.e2, freal);
fwrite(fid,xg.e3, freal);

fwrite(fid,xg.er, freal);    %spherical unit vectors
fwrite(fid,xg.etheta, freal);
fwrite(fid,xg.ephi, freal);

fwrite(fid,xg.r, freal);    %spherical coordinates
fwrite(fid,xg.theta, freal);
fwrite(fid,xg.phi, freal);

fwrite(fid,xg.x, freal);     %cartesian coordinates
fwrite(fid,xg.y, freal);
fwrite(fid,xg.z, freal);

fclose(fid);

end

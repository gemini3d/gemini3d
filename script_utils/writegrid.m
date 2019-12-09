function writegrid(xg, outdir, realbits)
%% write grid to raw binary files
% INCLUDes STUFF NOT NEEDED BY FORTRAN CODE BUT POSSIBLY USEFUL FOR PLOTTING

narginchk(2,3)
validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid parameters', 1)
validateattributes(outdir, {'char'}, {'vector'}, mfilename,'output directory',2)

if nargin < 3
  realbits = 64;
end
validateattr(realbits, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename,'real bits',3)

switch realbits
  case 64, freal = 'float64';
  case 32, freal = 'float32';
  otherwise, error(['unknown precision', num2str(realbits)])
end

%% MAKE THE OUTPUT DIRECTORY IF IT DOESN'T EXIST AND NOTIFY USER
outdir = absolute_path(outdir);
if ~isfolder(outdir)
  mkdir(outdir);
  disp(['Created: ', outdir])
end
% malformed paths can be "created" but are not accessible. Bug in Matlab mkdir().
assert(isfolder(outdir), [outdir, ' does not exist'])


filename = [outdir, '/simsize.dat'];
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

disp(['finishing writing grid to ', outdir])
end

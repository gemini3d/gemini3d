function writedata(ymd, UTsec,ns,vsx1,Ts, outdir, file_format, realbits)
%% WRITE STATE VARIABLE DATA TO BE USED AS INITIAL CONDITIONS
% FOR ANOTHER SIMULATION.  NOTE THAT WE
% DO NOT HERE OUTPUT ANY OF THE ELECTRODYNAMIC
% VARIABLES SINCE THEY ARE NOT NEEDED TO START THINGS
% UP IN THE FORTRAN CODE.
%
% INPUT ARRAYS SHOULD BE TRIMMED TO THE CORRECT SIZE
% I.E. THEY SHOULD NOT INCLUDE GHOST CELLS

narginchk(8,8)
validateattributes(ymd, {'numeric'}, {'vector', 'positive', 'numel', 3}, mfilename, 'year, month, day', 1)
validateattributes(UTsec, {'numeric'}, {'scalar', 'nonnegative'}, mfilename, 'seconds since UT midnight', 2)
validateattributes(ns, {'numeric'}, {'ndims', 4,'nonnegative'}, mfilename, 'density', 3)
validateattributes(vsx1, {'numeric'}, {'ndims', 4}, mfilename, 'velocity', 4)
validateattributes(Ts, {'numeric'}, {'ndims', 4,'nonnegative'}, mfilename, 'temperature', 5)
validateattributes(outdir, {'char'}, {'vector'}, mfilename, 'output directory',6)
validateattributes(file_format, {'char'}, {'vector'}, mfilename,'hdf5 or raw',7)
validateattributes(realbits, {'numeric'}, {'scalar','integer'},mfilename, '32 or 64',8)

outdir = absolute_path(outdir);
makedir(outdir)

switch file_format
  case {'h5','hdf5'}, write_hdf5(outdir, ymd, UTsec, ns, vsx1, Ts)
  case {'dat','raw'}, write_raw(outdir, ymd, UTsec, ns, vsx1, Ts, realbits)
  otherwise, error('unknown file_format')
end

end % function


function write_hdf5(outdir, ymd, UTsec, ns, vsx1, Ts)
fn = [outdir,'/initial_conditions.h5'];
disp(['write ',fn])
if isfile(fn), delete(fn), end

h5save(fn, '/ymd', int32(ymd))

freal = 'float32';

h5save(fn, '/UTsec', UTsec, [], freal)
h5save(fn, '/ns', ns, [], freal)
h5save(fn, '/vsx1', vsx1, [], freal)
h5save(fn, '/Ts', Ts, [], freal)

end % function


function write_raw(outdir, ymd, UTsec, ns, vsx1, Ts, realbits)

fn = [outdir,'/initial_conditions.dat'];
disp(['write ',fn])
fid=fopen(fn, 'w');

freal = ['float',int2str(realbits)];

fwrite(fid, ymd, freal);
fwrite(fid, UTsec, freal);
fwrite(fid,ns, freal);
fwrite(fid,vsx1, freal);
fwrite(fid,Ts, freal);

fclose(fid);

end

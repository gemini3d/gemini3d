function writedata(dmy,time,ns,vsx1,Ts,outdir,outID, format)

%--------------------------------------------------------
%-----WRITE STATE VARIABLE DATA TO BE USED AS INITIAL
%-----CONDITIONS FOR ANOTHER SIMULATION.  NOTE THAT WE
%-----DO NOT HERE OUTPUT ANY OF THE ELECTRODYNAMIC
%-----VARIABLES SINCE THEY ARE NOT NEEDED TO START THINGS
%-----UP IN THE FORTRAN CODE.
%-----
%-----INPUT ARRAYS SHOULD BE TRIMMED TO THE CORRECT SIZE
%-----(I.E. THEY SHOULD NOT INCLUDE GHOST CELLS
%--------------------------------------------------------
narginchk(8,8)
validateattributes(dmy, {'numeric'}, {'vector', 'positive', 'numel', 3}, mfilename, 'day, month, year', 1)
validateattributes(time, {'numeric'}, {'scalar', 'nonnegative'}, mfilename, 'time', 2)
validateattributes(ns, {'numeric'}, {'ndims', 4,'nonnegative'}, mfilename, 'density', 3)
validateattributes(vsx1, {'numeric'}, {'ndims', 4}, mfilename, 'velocity', 4)
validateattributes(Ts, {'numeric'}, {'ndims', 4,'nonnegative'}, mfilename, 'temperature', 5)
validateattributes(outdir, {'char'}, {'vector'}, mfilename, 'output directory',6)
validateattributes(outID, {'char'}, {'vector'}, mfilename, 'output filestem',7)
validateattr(format, {'char'}, {'vector'}, mfilename,'hdf5 or raw',8)

%% CHECK THAT THE DIRECTORY EXISTS - IF NOT CREATE IT AND GIVE THE USER AN INDICATION OF WHAT WAS DONE
if exist(outdir,'dir')~=7
  mkdir(outdir);
  disp(['Created: ', outdir])
end


switch format
  case 'hdf5', write_hdf5(outdir, outID, dmy, time, ns, vsx1, Ts)
  case 'raw', write_raw(outdir, outID, dmy, time, ns, vsx1, Ts)
end

end % function


function write_hdf5(outdir, outID, dmy, time, ns, vsx1, Ts)
fn = [outdir,'/',outID,'_ICs.h5'];
disp(['write ',fn])
h5save(fn, '/dmy', dmy)
h5save(fn, '/time', time)
h5save(fn, '/ns', ns)
h5save(fn, '/vsx1', vsx1)
h5save(fn, '/Ts', Ts)
end


function write_raw(outdir, outID, dmy, time, ns, vsx1, Ts)

fn = [outdir,'/',outID,'_ICs.dat'];
disp(['write ',fn])
fid=fopen(fn, 'w');
freal = 'float64';
fwrite(fid,dmy, freal);
fwrite(fid,time, freal);
fwrite(fid,ns, freal);
fwrite(fid,vsx1, freal);
fwrite(fid,Ts, freal);

fclose(fid);

end

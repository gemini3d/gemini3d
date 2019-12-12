function writedata(dmy,time,ns,vsx1,Ts,outdir,outID, realbits)

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
narginchk(7,8)
validateattributes(dmy, {'numeric'}, {'vector', 'positive', 'numel', 3}, mfilename, 'day, month, year', 1)
validateattributes(time, {'numeric'}, {'scalar', 'nonnegative'}, mfilename, 'time', 2)
validateattributes(ns, {'numeric'}, {'ndims', 4,'nonnegative'}, mfilename, 'density', 3)
validateattributes(vsx1, {'numeric'}, {'ndims', 4}, mfilename, 'velocity', 4)
validateattributes(Ts, {'numeric'}, {'ndims', 4,'nonnegative'}, mfilename, 'temperature', 5)
validateattributes(outdir, {'char'}, {'vector'}, mfilename, 'output directory',6)
validateattributes(outID, {'char'}, {'vector'}, mfilename, 'output filestem',7)

if nargin < 8
  realbits = 64;
end
validateattr(realbits, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename,'real bits',8)

switch realbits
  case 64, freal = 'float64';
  case 32, freal = 'float32';
  otherwise, error(['unknown precision', num2str(realbits)])
end

%% CHECK THAT THE DIRECTORY EXISTS - IF NOT CREATE IT AND GIVE THE USER AN INDICATION OF WHAT WAS DONE
if exist(outdir,'dir')~=7
  mkdir(outdir);
  disp(['Created: ', outdir])
end

%% OPEN THE FILE
filename=[outdir,filesep,outID,'_ICs.dat'];
fid=fopen(filename,'w');

%% WRITE THE DATA
fwrite(fid,dmy, freal);
fwrite(fid,time, freal);
fwrite(fid,ns, freal);
fwrite(fid,vsx1, freal);
fwrite(fid,Ts, freal);

fclose(fid);

end % function

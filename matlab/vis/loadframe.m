function dat = loadframe(direc,ymd,UTsec,flagoutput,mloc,xg,file_format, config_file, realbits)

narginchk(3,9)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)
validateattr(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 2)
validateattr(UTsec, {'numeric'}, {'vector'}, mfilename, 'UTC second', 3)

if nargin < 9 || isempty(realbits), realbits = 64; end

if nargin < 8 || isempty(config_file)
  config_file = [direc, '/inputs'];
end

if nargin < 5 || isempty(flagoutput) || isempty(mloc)
  p = read_config(config_file);
  flagoutput = p.flagoutput;
  mloc = p.mloc;
end
validateattr(flagoutput,{'numeric'},{'scalar'},mfilename,'output flag',4)

if ~isempty(mloc)
  validateattr(mloc, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'magnetic coordinates', 5)
end

if nargin < 7 || isempty(file_format)
  file_format = 'auto';
end
validateattr(file_format, {'char'}, {'vector'}, mfilename, 'raw or hdf5', 7)

if nargin < 6 || isempty(xg)
  xg = readgrid([direc, '/inputs'], file_format, realbits);
end
validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 6)


direc = absolute_path(direc);

%% LOAD DIST. FILE
stem0 = datelab(ymd, UTsec);
switch file_format
  case {'h5','hdf5'}, suffix = {'.h5'};
  case {'dat','raw'}, suffix = {'.dat'};
  case {'nc'}, suffix = {'.nc'};
  otherwise, suffix = {'.h5', '.nc', '.dat'};
end

for ext = suffix
  stem = stem0;
  filename = [direc, filesep, stem, ext{1}];
  if is_file(filename)
    break
  end
  % switch microsecond to one for first time step
  stem(end) = '1';
  filename = [direc, filesep,stem, ext{1}];
  if is_file(filename)
    break
  end
end

assert(is_file(filename), ['loadframe: ', filename, ' does not exist.'])

switch flagoutput
  case 1
    dat = loadframe3Dcurv(filename);
  case 2
    dat = loadframe3Dcurvavg(filename);
  otherwise
    dat = loadframe3Dcurvne(filename);
end

%% SET MAGNETIC LATITUDE AND LONGITUDE OF THE SOURCE
dat.mlatsrc=[];
dat.mlonsrc=[];
if nargout >= 2 && ~isempty(mloc)
  dat.mlatsrc=mloc(1);
  dat.mlonsrc=mloc(2);
end

end % function

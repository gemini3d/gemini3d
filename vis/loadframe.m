function [ne,mlatsrc,mlonsrc,xg,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = loadframe(direc,ymd,UTsec,flagoutput,mloc,xg,file_format, config_file)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, '/../script_utils'])

narginchk(3,8)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)
validateattr(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 2)
validateattr(UTsec, {'numeric'}, {'vector'}, mfilename, 'UTC second', 3)

if nargin < 8 || isempty(config_file)
  config_file = [direc, '/inputs'];
end

if nargin < 5 || isempty(flagoutput) || isempty(mloc)
  [~,~,~,~,flagoutput,mloc] = readconfig(config_file);
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
  xg = readgrid([direc, '/inputs'], file_format);
end
validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 6)


direc = absolute_path(direc);

%% SET MAGNETIC LATITUDE AND LONGITUDE OF THE SOURCE
if nargout >= 2 && ~isempty(mloc)
  mlatsrc=mloc(1);
  mlonsrc=mloc(2);
else
  mlatsrc=[];
  mlonsrc=[];
end

%% LOAD DIST. FILE
stem0 = datelab(ymd, UTsec);
switch file_format
  case 'hdf5', suffix = {'.h5'};
  case 'raw', suffix = {'.dat'};
  otherwise, suffix = {'.h5', '.dat'};
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
    [ne,v1,Ti,Te,J1,v2,v3,J2,J3,ns,vs1,Ts,Phitop] = loadframe3Dcurv(filename);
  case 2
    [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg(filename);
    ns=[]; vs1=[]; Ts=[];
  otherwise
    ne = loadframe3Dcurvne(filename);
    v1=[]; Ti=[]; Te=[]; J1=[]; v2=[]; v3=[]; J2=[]; J3=[];
    ns=[]; vs1=[]; Ts=[]; Phitop=[];
end

end % function

function [ne,mlatsrc,mlonsrc,xg,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = loadframe(direc,ymd,UTsec,flagoutput,mloc,xg)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, '..', filesep, 'script_utils'])

narginchk(3,6)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)
validateattr(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 2)
validateattr(UTsec, {'numeric'}, {'vector'}, mfilename, 'UTC second', 3)

if nargin < 5 || isempty(flagoutput) || isempty(mloc)
  [~,~,~,~,flagoutput,mloc] = readconfig([direc, filesep, 'inputs']);
end
validateattr(flagoutput,{'numeric'},{'scalar'},mfilename,'output flag',4)
validateattr(mloc, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'magnetic coordinates', 5)

if nargin < 6 || isempty(xg)
  xg = readgrid([direc, filesep, 'inputs']);
end
validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 6)

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
for ext = {'.h5', '.dat'}
  stem = stem0;
  filename = [stem, ext{1}];
  if is_file([direc, filesep, filename])
    break
  end
  % switch microsecond to one for first time step
  stem(end) = '1';
  filename = [stem, ext{1}];
  if is_file([direc, filesep, filename])
    break
  end
end

assert(is_file([direc, filesep, filename]), ['loadframe: ', [direc,filesep,stem,'.{dat,h5}'], ' does not exist.'])

switch flagoutput
  case 1
    [ne,v1,Ti,Te,J1,v2,v3,J2,J3,ns,vs1,Ts,Phitop] = loadframe3Dcurv(direc,filename);
  case 2
    [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg(direc,filename);
    ns=[]; vs1=[]; Ts=[];
  otherwise
    ne=loadframe3Dcurvne(direc,filename);
    v1=[]; Ti=[]; Te=[]; J1=[]; v2=[]; v3=[]; J2=[]; J3=[];
    ns=[]; vs1=[]; Ts=[]; Phitop=[];
end

end % function

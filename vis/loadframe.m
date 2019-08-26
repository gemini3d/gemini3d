function [ne,mlatsrc,mlonsrc,xg,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = loadframe(direc,ymd,UTsec,ymd0,UTsec0,tdur,dtout,flagoutput,mloc,xg)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, '..', filesep, 'script_utils'])

narginchk(3,10)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)
validateattr(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 2)
validateattr(UTsec, {'numeric'}, {'vector'}, mfilename, 'UTC second', 3)

if nargin>=4
  validateattr(ymd0, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 4)
end
if nargin>=5
  validateattr(UTsec0, {'numeric'}, {'scalar'}, mfilename, 'UTC second', 5)
end

if nargin>=6 && ~isempty(tdur)
  validateattr(tdur,{'numeric'},{'scalar'},mfilename,'simulation duration',6)
end
if nargin>=7 && ~isempty(dtout)
  validateattr(dtout,{'numeric'},{'scalar'},mfilename,'output time step',7)
end

if nargin>=8 && ~isempty(flagoutput)
  validateattr(flagoutput,{'numeric'},{'scalar'},mfilename,'output flag',8)
end
if nargin>=9 && ~isempty(mloc)
  validateattr(mloc, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'magnetic coordinates', 9)
end
if nargin>=10 && ~isempty(xg)
  validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 10)
end


% READ IN THE SIMULATION INFORMATION IF IT HAS NOT ALREADY BEEN PROVIDED
if (~exist('UTsec0','var') || ~exist('ymd0','var') || ~exist('mloc','var') || ~exist('tdur','var') ...
     || ~exist('dtout','var') || ~exist('flagoutput','var') || ~exist('mloc','var') )
  [ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs/config.ini']);
end


% CHECK WHETHER WE NEED TO RELOAD THE GRID (WHICH CAN BE TIME CONSUMING)
if nargout >= 4 && ~exist('xg','var')
  xg = readgrid([direc, filesep, 'inputs']);
end


%% SET MAGNETIC LATITUDE AND LONGITUDE OF THE SOURCE
if nargout >= 2 && ~isempty(mloc)
  mlatsrc=mloc(1);
  mlonsrc=mloc(2);
else
  mlatsrc=[];
  mlonsrc=[];
end


%% LOAD DIST. FILE
stem = datelab(ymd, UTsec);
filename = [stem, '.dat'];
if ~exist(filename, 'file') % switch microsecond to one for first time step
    stem(end) = '1';
    filename = [stem, '.dat'];
end

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

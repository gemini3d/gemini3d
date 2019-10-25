function [ne,mlatsrc,mlonsrc,xg,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = loadframe(direc,ymd,UTsec,flagoutput,mloc,xg)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, '..', filesep, 'script_utils'])

narginchk(3,6)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)
validateattr(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 2)
validateattr(UTsec, {'numeric'}, {'vector'}, mfilename, 'UTC second', 3)

if nargin>=4 && ~isempty(flagoutput)
  validateattr(flagoutput,{'numeric'},{'scalar'},mfilename,'output flag',4)
end
if nargin>=5 && ~isempty(mloc)
  validateattr(mloc, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'magnetic coordinates', 5)
end
if nargin>=6 && ~isempty(xg)
  validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 6)
end


% READ IN THE SIMULATION INFORMATION IF IT HAS NOT ALREADY BEEN PROVIDED
if ~exist('flagoutput','var') || ~exist('mloc','var')
  [~,~,~,~,flagoutput,mloc] = readconfig([direc, filesep, 'inputs']);
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
if ~is_file([direc, filesep, filename]) % switch microsecond to one for first time step
    filenameold = filename;
    stem(end) = '1';
    filename = [stem, '.dat'];
    assert(is_file([direc, filesep, filename]), ['loadframe: ', [direc,filesep,filenameold], ' does not exist.'])
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

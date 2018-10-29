function [ne,mlatsrc,mlonsrc,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = loadframe(direc,UTsec,ymd,UTsec0,ymd0,mloc,xg)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd,'/../script_utils'])

narginchk(5, 7)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)
validateattr(UTsec, {'numeric'}, {'vector'}, mfilename, 'UTC second', 2)
validateattr(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 3)
validateattr(UTsec, {'numeric'}, {'scalar'}, mfilename, 'UTC second', 4)
validateattr(ymd0, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 5)
if nargin>=6 && ~isempty(mloc)
  validateattr(mloc, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'magnetic coordinates', 6)
end
%if nargin<7
%  xg = [];
%else
%  validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 7)
%end
%% READ IN THE SIMULATION INFORMATION IF IT HAS NOT ALREADY BEEN PROVIDED
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs/config.ini']);

%% CHECK WHETHER WE NEED TO RELOAD THE GRID (WHICH CAN BE TIME CONSUMING)
%if isempty(xg)
%  xg = readgrid([direc,'/inputs/']);
%end

%% SET MAGNETIC LATITUDE AND LONGITUDE OF THE SOURCE
if ~isempty(mloc)
  mlatsrc=mloc(1);
  mlonsrc=mloc(2);
else
  mlatsrc=[];
  mlonsrc=[];
end


%% LOAD DIST. FILE
filestr=datelab(ymd,UTsec);
if ymd(1)==ymd0(1) && ymd(2)==ymd0(2) && ymd(3)==ymd0(3) && UTsec==UTsec0    %tack on the decimal part
%  filename=[filestr,'.000001.dat']
  filestr(end)='1';
%else
%  filename=[filestr,'.000000.dat']
end
filename=[filestr,'.dat'];

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

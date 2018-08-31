function [ne,mlatsrc,mlonsrc,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = loadframe(direc,UTsec,ymd,UTsec0,ymd0,autoload,flagoutput,mloc,xg)

%PATH TO PLOTTING FUNCTIONS AND SHARED SCRIPT UTILITIES
addpath('../script_utils')


%READ IN THE SIMULATION INFORMATION IF IT HAS NOT ALREADY BEEN PROVIDED
if (~exist('ymd0','var'))
  [ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs/config.ini']);
end 


%SET GRID AUTOLOADING - needs to be a function input
if(~exist('autoload','var'))
  autoload=1;
end


%CHECK WHETHER WE NEED TO RELOAD THE GRID (WHICH CAN BE TIME CONSUMING)
if (~exist('xg','var') & autoload)
  %WE ALSO NEED TO LOAD THE GRID FILE
  xg=readgrid([direc,'/inputs/']);
end


%SET THE MAGNETIC LATITUDE AND LONGITUDE OF THE SOURCE
if (~isempty(mloc))
  mlatsrc=mloc(1);
  mlonsrc=mloc(2);
else
  mlatsrc=[];
  mlonsrc=[];
end


%LOAD DIST. FILE
filestr=datelab(ymd,UTsec);
if (ymd(1)==ymd0(1) & ymd(2)==ymd0(2) & ymd(3)==ymd0(3) & UTsec==UTsec0)      %tack on the decimal part
%  filename=[filestr,'.000001.dat']
  filestr(end)='1';
%else
%  filename=[filestr,'.000000.dat']
end
filename=[filestr,'.dat'];

if (flagoutput==1)
  [ne,v1,Ti,Te,J1,v2,v3,J2,J3,ns,vs1,Ts,Phitop] = loadframe3Dcurv(direc,filename);
elseif (flagoutput==2)
  [ne,v1,Ti,Te,J1,v2,v3,J2,J3,Phitop] = loadframe3Dcurvavg(direc,filename);
  ns=[]; vs1=[]; Ts=[];
else
   ne=loadframe3Dcurvne(direc,filename);
   v1=[]; Ti=[]; Te=[]; J1=[]; v2=[]; v3=[]; J2=[]; J3=[];
   ns=[]; vs1=[]; Ts=[]; Phitop=[];
end

end % function

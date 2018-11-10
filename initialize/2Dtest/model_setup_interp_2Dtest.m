%LOWRES 2D EXAMPLE FOR TESTING
xdist=200e3;    %eastward distance
ydist=600e3;    %northward distance
lxp=80;
lyp=1;
glat=67.11;
glon=212.95;
gridflag=0;
I=90;


%ADD PATHS FOR FUNCTIONS
cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, '..', filesep,'..',filesep,'script_utils']);
addpath([cwd, filesep, '..', filesep,'..',filesep,'setup']);
addpath([cwd, filesep, '..', filesep,'..',filesep,'setup',filesep,'gridgen'])
addpath([cwd, filesep, '..', filesep,'..',filesep,'vis']);


%RUN THE GRID GENERATION CODE
if (~exist('xg'))
  xg=makegrid_cart_3D(xdist,lxp,ydist,lyp,I,glat,glon);
end

eqdir='../../../simulations/2Dtest_eq/';
distdir='../../../simulations/2Dtest/';
simID='2Dtest';
[nsi,vs1i,Tsi,xgin,ns,vs1,Ts]=eq2dist(eqdir,distdir,simID,xg);

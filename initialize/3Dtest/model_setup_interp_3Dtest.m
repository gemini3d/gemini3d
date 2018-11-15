%LOWRES 2D EXAMPLE FOR TESTING
xdist=200e3;    %eastward distance
ydist=200e3;    %northward distance
lxp=20;
lyp=20;
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
  xg=makegrid_cart_3D_lowresx1(xdist,lxp,ydist,lyp,I,glat,glon);
end


%PRODUCE AND IC FROM EQ SIMULATION
eqdir='../../../simulations/3Dtest_eq/';
simID='3Dtest';
[nsi,vs1i,Tsi,xgin,ns,vs1,Ts]=eq2dist(eqdir,simID,xg);

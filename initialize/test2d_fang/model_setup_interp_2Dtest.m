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
gemdir = [cwd, '/../..'];
assert(isfolder(gemdir), [gemdir, ' not found'])
addpath([gemdir, '/script_utils']);
addpath([gemdir, '/setup']);
addpath([gemdir, '/setup/gridgen'])
addpath([gemdir, '/vis']);


%RUN THE GRID GENERATION CODE
xg=makegrid_cart_3D(xdist,lxp,ydist,lyp,I,glat,glon);

eqdir = '../../../simulations/2Dtest_eq';
simID='2Dtest';
[nsi,vs1i,Tsi,xgin,ns,vs1,Ts]=eq2dist(eqdir,simID,xg);

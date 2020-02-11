 %% LOWRES 2D EXAMPLE FOR TESTING
 p.format = 'hdf5';
 p.xdist=1200e3;    %eastward distance
 p.ydist=600e3;    %northward distance
 p.lxp=15;
 p.lyp=15;
 p.glat=67.11;
 p.glon=212.95;
 p.I=90;


%ADD PATHS FOR FUNCTIONS
cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, '..', filesep,'..',filesep,'setup']);
addpath([cwd, filesep, '..', filesep,'..',filesep,'setup',filesep,'gridgen'])

%% GRID GENERATION
xg = makegrid_cart_3D_lowresx1(p);


%GENERATE SOME INITIAL CONDITIONS FOR A PARTICULAR EVENT
%%ISINGLASS B LAUNCH
UT=7.5;
dmy=[2,3,2017];
activ=[76.5,79.3,31.5];


%USE OLD CODE FROM MATLAB MODEL
nmf=5e11;
nme=2e11;
[ns,Ts,vsx1]=eqICs3D(xg,UT,dmy,activ,nmf,nme);    %note that this actually calls msis_matlab - should be rewritten to include the neutral module form the fortran code!!!


%WRITE THE GRID AND INITIAL CONDITIONS
outdir = [cwd, '/../../gemini_sim/test3d_eq'];
writegrid(xg,outdir);
time=UT*3600;   %doesn't matter for input files
writedata(dmy,time,ns,vsx1,Ts,outdir,p.format);

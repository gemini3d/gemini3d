%RISR LOWRES GRID (CARTESIAN)
xdist=1.5e6;
ydist=1.5e6;
lxp=20;
lyp=20;
glat=75.6975;
glon=360.0-94.8322;
gridflag=0;
I=90;


%ADD PATHS TO SCRIPT UTILS (hopefully this catches the fact that the makegrids need access to geomag conversion functions
addpath ../../script_utils;
addpath ../../setup/gridgen;
addpath ../../setup;

%MATLAB GRID GENERATION
xg=makegrid_cart_3D_lowresx1(xdist,lxp,ydist,lyp,I,glat,glon);


%GENERATE SOME INITIAL CONDITIONS FOR A PARTICULAR EVENT
%RISR
UT=18000/3600;
dmy=[20,2,2013];
activ=[150.0,150.0,50.0];


%USE OLD CODE FROM MATLAB MODEL
nmf=5e11;
nme=2e11;
[ns,Ts,vsx1]=eqICs3D(xg,UT,dmy,activ,nmf,nme);    %note that this actually calls msis_matlab - should be rewritten to include the neutral module form the fortran code!!!


%WRITE THE GRID AND INITIAL CONDITIONS
outdir='~/zettergmdata/simulations/input/RISR_eq/'
simlabel='RISR_eq'
writegrid(xg,outdir);
time=UT*3600;   %doesn't matter for input files
writedata(dmy,time,ns,vsx1,Ts,outdir,simlabel);


%RESET PATH
rmpath ../../script_utils;
rmpath ../../setup/gridgen;
rmpath ../../setup;

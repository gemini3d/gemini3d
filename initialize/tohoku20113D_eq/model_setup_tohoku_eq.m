%A HIGHRES TOHOKU INIT GRID
dtheta=10;
dphi=15;
lp=100;
lq=500;
lphi=25;
altmin=80e3;
glat=42.45;
glon=143.4;
gridflag=1;
flagsource=1;


%ADD PATHS TO SCRIPT UTILS (hopefully this catches the fact that the makegrids need access to geomag conversion functions
addpath ../../script_utils;
addpath ../../setup/gridgen;
addpath ../../setup;

%MATLAB GRID GENERATION
xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);


%GENERATE SOME INITIAL CONDITIONS FOR A PARTICULAR EVENT, THESE ACTUALLY DON'T MATTER MUCH SO YOU CAN MAKE UP STUFF
UT=5.75;
dmy=[11,3,2011];
activ=[120,120,25];


%USE OLD CODE FROM MATLAB MODEL
nmf=5e11;
nme=2e11;
[ns,Ts,vsx1]=eqICs3D(xg,UT,dmy,activ,nmf,nme);    %note that this actually calls msis_matlab - should be rewritten to include the neutral module form the fortran code!!!


%WRITE THE GRID AND INITIAL CONDITIONS
outdir='~/'
simlabel='tohoku_eq'
writegrid(xg,outdir,simlabel);
time=UT*3600;   %doesn't matter for input files
writedata(dmy,time,ns,vsx1,Ts,outdir,simlabel);


%RESET PATH
rmpath ../../script_utils;
rmpath ../../setup/gridgen;
rmpath ../../setup;


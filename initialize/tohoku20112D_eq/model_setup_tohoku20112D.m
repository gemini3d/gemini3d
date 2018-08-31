%EQ TOHOKU GRID 2D
dtheta=10;
dphi=12;
glat=42.45;
glon=143.40;
lp=256;
lq=750;
lphi=1;
altmin=80e3;
gridflag=1;


%ADD PATHS TO SCRIPT UTILS (hopefully this catches the fact that the makegrids need access to geomag conversion functions
addpath ../../script_utils;
addpath ../../setup/gridgen;
addpath ../../setup;

%MATLAB GRID GENERATION
xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);


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
outdir='~/zettergmdata/simulations/input/tohoku20112D_eq/'
simlabel='tohoku20112D_eq'
writegrid(xg,outdir);
time=UT*3600;   %doesn't matter for input files
writedata(dmy,time,ns,vsx1,Ts,outdir,simlabel);


%RESET PATH
rmpath ../../script_utils;
rmpath ../../setup/gridgen;
rmpath ../../setup;

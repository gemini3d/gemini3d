 %NEPAL 2015 GRID
 dtheta=10;
 dphi=16;
 lp=256;
 lq=750;
 lphi=1;
 altmin=80e3;
 glat=35.75;
 glon=84.73;
 gridflag=1;


%ADD PATHS TO SCRIPT UTILS (hopefully this catches the fact that the makegrids need access to geomag conversion functions
addpath ../../script_utils;
addpath ../../setup/gridgen;
addpath ../../setup;

%MATLAB GRID GENERATION
xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);


%GENERATE SOME INITIAL CONDITIONS FOR A PARTICULAR EVENT
%NEPAL, 2015
activ=[136,125.6,0.5];
dmy=[25,4,2015];
t0=(6+11/60)*3600;
UT=t0/3600;


%USE OLD CODE FROM MATLAB MODEL
nmf=5e11;
nme=2e11;
[ns,Ts,vsx1]=eqICs3D(xg,UT,dmy,activ,nmf,nme);    %note that this actually calls msis_matlab - should be rewritten to include the neutral module form the fortran code!!!


%WRITE THE GRID AND INITIAL CONDITIONS
simlabel='nepal20152D_eq'
outdir='~/zettergmdata/simulations/input/nepal20152D_eq'
writegrid(xg,outdir);
time=UT*3600;   %doesn't matter for input files
writedata(dmy,time,ns,vsx1,Ts,outdir,simlabel);


%RESET PATH
rmpath ../../script_utils;
rmpath ../../setup/gridgen;
rmpath ../../setup;

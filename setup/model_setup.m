% %HIGH-LATITUDE GRID (FOR TESTING ELECTRODYNAMICS)
% dtheta=2;
% dphi=5;
% lp=100;
% lq=200;
% lphi=100;
% altmin=80e3;
% glat=65;    %high-latitude
% glon=270;
% gridflag=0;

%HIGH-LATITUDE GRID (COARSE)
% dtheta=2;
% dphi=5;
% lp=35;
% lq=200;
% lphi=44;
% altmin=80e3;
% glat=65;    %high-latitude
% glon=270;
% gridflag=0;

% %LOW-LAT. GRID (COARSE)
% dtheta=10;
% dphi=10;
% lp=35;
% lq=350;
% lphi=48;
% altmin=80e3;
% glat=25;    %low-latitude
% glon=270;
% gridflag=1;

%LOW-LAT. GRID (Higher res...)
%dtheta=10;
%dphi=10;
%lp=250;
%lq=1000;
%lphi=240;
%altmin=80e3;
%glat=25;    %low-latitude
%glon=270;
%gridflag=1;


%3D TOHOKU GRID (VERY COARSE)
%dtheta=7.5;
%dphi=7.5;
%glat=42.45;
%glon=143.40;
%lp=150;
%lq=500;
%lphi=150;
%altmin=80e3;
%gridflag=1;


%EQ TOHOKU GRID
%dtheta=7.5;
%dphi=12;
%glat=42.45;
%glon=143.40;
%lp=75;
%lq=500;
%lphi=15;
%altmin=80e3;
%gridflag=1;


%%EQ TOHOKU GRID 2D
%dtheta=7.5;
%dphi=12;
%glat=42.45;
%glon=143.40;
%lp=256;
%lq=750;
%lphi=1;
%altmin=80e3;
%gridflag=1;


%%EQ TOHOKU GRID 2D
%dtheta=7.5;
%dphi=12;
%glat=42.45;
%glon=142.73;
%lp=256;
%lq=750;
%lphi=1;
%altmin=80e3;
%gridflag=1;
%

%
%%EQ TOHOKU GRID 2D WIDE
%dtheta=11;
%dphi=12;
%glat=42.45;
%glon=142.73;
%lp=20*16;
%lq=750;
%lphi=1;
%altmin=80e3;
%gridflag=1;

%  %CHILE 2015 GRID
%  dtheta=9;
%  dphi=14;
%  lp=350;
%  lq=500;
%  lphi=1;
%  altmin=80e3;
%  glat=17.0;
%  glon=288.2;
%  gridflag=1;
%
% %NEPAL 2015 GRID
% dtheta=10;
% dphi=16;
% lp=350;
% lq=500;
% lphi=1;
% altmin=80e3;
% glat=35.75;
% glon=84.73;
% gridflag=1;
%
% %MOORE, OK GRID (FULL)
% dtheta=25;
% dphi=35;
%% lp=800;
% lp=256;
%% lq=1250;
% lq=400;
% lphi=1;
% altmin=80e3;
% glat=39;
% glon=262.51;
% gridflag=0;
%

%%MOORE, OK GRID (PARTIAL)
%dtheta=15;
%dphi=20;
%lp=100;
%lq=350;
%lphi=100;
%altmin=80e3;
%glat=39;
%glon=262.51;
%gridflag=0;

%%RISR LOWRES GRID (CARTESIAN)
%xdist=1.5e6;
%ydist=1.5e6;
%lxp=10;
%lyp=10;
%glat=75.6975;
%glon=360.0-94.8322;
%gridflag=0;
%I=90;
%

%{
%CHILE 2015 GRID
dtheta=9;
dphi=14;
lp=350;
lq=500;
lphi=1;
altmin=80e3;
glat=19.25;
glon=287.5;
gridflag=1;
%}

%{
%PFISR LOWRES GRID (CARTESIAN)
xdist=640e3;    %eastward distance
ydist=385e3;    %northward distance
lxp=150;
lyp=192;
glat=67.11;
glon=212.95;
gridflag=0;
I=90;
%}
%
%%PFISR LOWRES GRID (CARTESIAN)
%xdist=1200e3;    %eastward distance
%ydist=600e3;    %northward distance
%lxp=15;
%lyp=15;
%glat=67.11;
%glon=212.95;
%gridflag=0;
%I=90;
%

%%LOWRES 2D EXAMPLE FOR TESTING
%xdist=1200e3;    %eastward distance
%ydist=600e3;    %northward distance
%lxp=15;
%lyp=1;
%glat=67.11;
%glon=212.95;
%gridflag=0;
%I=90;
%


% %LOWRES 2D EXAMPLE FOR TESTING
% xdist=1200e3;    %eastward distance
% ydist=600e3;    %northward distance
% lxp=1;
% lyp=15;
% glat=67.11;
% glon=212.95;
% gridflag=0;
% I=90;


% %EXAMPLE FOR KRISTINA MIDEX MISSION
% xdist=700e3;    %eastward distance
% ydist=150e3;    %northward distance
% lxp=20;
% lyp=20;
% glat=67.11;
% glon=212.95;
% gridflag=0;
% I=90;


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
addpath ../script_utils;
addpath ./gridgen;


%MATLAB GRID GENERATION
xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
%xg=makegrid_cart_3D(xdist,lxp,ydist,lyp,I,glat,glon);


%GENERATE SOME INITIAL CONDITIONS FOR A PARTICULAR EVENT

%NEPAL, 2015
%activ=[136,125.6,0.5];
%dmy=[25,4,2015];
%t0=(6+11/60)*3600;
%UT=t0/3600;

%RISR
%UT=18000/3600;
%dmy=[20,2,2013];
%activ=[150.0,150.0,50.0];

%CHILE 2015
%UT=82473;
%dmy=[16,9,2015];
%activ=[109,109,5];

%%ISINGLASS B LAUNCH
UT=7.5;
dmy=[2,3,2017];
activ=[76.5,79.3,31.5];


%USE OLD CODE FROM MATLAB MODEL
nmf=5e11;
nme=2e11;
[ns,Ts,vsx1]=eqICs3D(xg,UT,dmy,activ,nmf,nme);    %note that this actually calls msis_matlab - should be rewritten to include the neutral module form the fortran code!!!


%WRITE THE GRID AND INITIAL CONDITIONS
%simlabel='chile2015_eq'
%simlabel='3DPCarc_eq'
%simlabel='mooreOKfull_eq'
%simlabel='nepal2D_eq'
%simlabel='2Dtest_nonperumuted_eq'
%simlabel='ARCS_eq'
outdir='~/'
simlabel='tohoku_eq'
writegrid(xg,outdir,simlabel);
time=UT*3600;   %doesn't matter for input files
writedata(dmy,time,ns,vsx1,Ts,outdir,simlabel);


%RESET PATH
rmpath ../script_utils;
rmpath ./gridgen;


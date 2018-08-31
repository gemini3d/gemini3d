%SUPER LOWRES
%{
dtheta=7.5;
dphi=11;
glat=42.45;
glon=143.40;
lp=75;
%lp=150;
lq=500;
%lq=750;
lphi=80;
altmin=80e3;
gridflag=1;
%}


%{
%SOMETHING HIGHER RES FOR CELERITY NODES
dtheta=7.25;
dphi=10.75;
glat=42.45;
glon=143.40;
lp=180;
lq=1000;
lphi=180;
altmin=80e3;
gridflag=1;



%MATLAB GRID GENERATION
if (~exist('xg'))
  xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
end

%{
%GENERATE SOME INITIAL CONDITIONS
UT=0;
dmy=[11,3,2011];
activ=[100,100,10];
nmf=5e11;
nme=2e11;
[ns,Ts,vsx1]=eqICs3D(xg,UT,dmy,activ,nmf,nme);    %note that this actually calls msis_matlab - should be rewritten to include the neutral module form the fortran code!!!
%}
%}


%%UNIFORM HI-RES EXAMPLE
%dtheta=7.25;
%dphi=10.75;
%lp=288;
%lq=500;
%lphi=288;
%altmin=80e3;
%glat=42.45;
%glon=143.4;
%gridflag=1;
%

%%NONUNIFORM GRID
%%THIS WILL ONLY WORK ON LEVITY (NEEDS MORE THAN 64 GB RAM)
%dtheta=7.05;
%dphi=10.5;
%lp=320;
%lq=50;
%lphi=320;
%altmin=80e3;
%glat=42.45;
%glon=143.4;
%gridflag=1;


%%2D EXAMPLE, HIGHRES
%dtheta=7.25;
%dphi=12;
%glat=42.45;
%glon=142.73;
%lp=1024;
%lq=1024;
%lphi=1;
%altmin=80e3;
%gridflag=1;

%  %CHILE 2015 GRID
%  dtheta=8;
%  dphi=14;
% % lp=1024;
% % lq=1024;
%  lp=350;
%  lq=350;
%  lphi=1;
%  altmin=80e3;
%  glat=17.0;
%  glon=288.2;
%  gridflag=1;

%%NEPAL 2015 GRID
%dtheta=8;
%dphi=14;
%lp=1024;
%lq=1024;
%lphi=1;
%altmin=80e3;
%glat=35.75;
%glon=84.73;
%gridflag=1;
%
%
% %MOORE, OK GRID (FULL)
% dtheta=24.5;
% dphi=34.5;
% lp=800;
%% lp=256;
% lq=1250;
%% lq=400;
% lphi=1;
% altmin=80e3;
% glat=39;
% glon=262.51;
% gridflag=0;


%%RISR LOWRES GRID (CARTESIAN)
%xdist=1e6;
%ydist=1e6;
%lxp=125;
%lyp=200;
%glat=75.6975;
%glon=360.0-94.8322;
%gridflag=0;
%I=90;
%

%
%%RISR PERIODIC GDI RUN
%xdist=200e3;
%ydist=200e3;
%lxp=250;
%lyp=250;
%glat=75.6975;
%glon=360.0-94.8322;
%gridflag=0;
%I=90;
%

%
%%RISR PERIODIC GDI RUN (HIGHRES)
%xdist=200e3;
%ydist=200e3;
%lxp=1000;
%lyp=1000;
%glat=75.6975;
%glon=360.0-94.8322;
%gridflag=0;
%I=90;


%%RISR PERIODIC KHI RUN
%xdist=40e3;
%ydist=22e3;
%lxp=200;
%lyp=110;
%glat=75.6975;
%glon=360.0-94.8322;
%gridflag=0;
%I=90;
%

%{
%RISR PERIODIC KHI RUN
xdist=40e3;
ydist=22e3;
lxp=400;
lyp=220;
glat=75.6975;
glon=360.0-94.8322;
gridflag=0;
I=90;
%}

% %CHILE 2015 GRID 3D
% dtheta=8.5;
% dphi=12;
% lp=200;
% lq=450;
% lphi=128;
% altmin=80e3;
% glat=17.0;
% glon=288.2;
% gridflag=1;
%
%
%%PFISR LOWRES GRID (CARTESIAN)
%xdist=500e3;    %eastward distance
%ydist=250e3;    %northward distance
%lxp=100;
%lyp=200;
%glat=67.11;
%glon=212.95;
%gridflag=0;
%I=90;
%

% %CHILE 2015 GRID 3D (highres)
%  dtheta=8.5;
%  dphi=12;
%  lp=440;
%  lq=550;
%  lphi=440;
%  altmin=80e3;
%  glat=17.0;
%  glon=288.2;
%  gridflag=1;
% 
% 
% %CHILE 2015 GRID 3D (veryhighres)
%  dtheta=8.5;
%  dphi=12;
%  lp=550;
%  lq=650;
%  lphi=512;
%  altmin=80e3;
%  glat=17.0;
%  glon=288.2;
%  gridflag=1;

%{
%CHILE 2015 GRID 3D (medres)
 dtheta=8.5;
 dphi=12;
 lp=440;
 lq=550;
 lphi=352;
 altmin=80e3;
 glat=17.0;
 glon=288.2;
 gridflag=1;
%}


%{
%LOWRES 2D EXAMPLE FOR TESTING
xdist=200e3;    %eastward distance
ydist=600e3;    %northward distance
lxp=80;
lyp=1;
glat=67.11;
glon=212.95;
gridflag=0;
I=90;
%}

% %PFISR LOWRES GRID (CARTESIAN)
% xdist=500e3;    %eastward distance
% ydist=100e3;    %northward distance
% lxp=384/2;
% lyp=512/2;
% glat=67.11;
% glon=212.95;
% gridflag=0;
% I=90;

%{
%A HIGHRES TOHOKU
dtheta=7.5;
dphi=12;
lp=576;
lq=580;
lphi=576;
altmin=80e3;
glat=42.45;
glon=143.4;
gridflag=1;
flagsource=1;
%}


%A MEDIUM RES TOHOKU
dtheta=7.5;
dphi=12;
lp=256;
lq=580;
lphi=288;
altmin=80e3;
glat=42.45;
glon=143.4;
gridflag=1;
flagsource=1;


%{
%A lowRES TOHOKU
dtheta=7.5;
dphi=12;
lp=128;
lq=580;
lphi=150;
altmin=80e3;
glat=42.45;
glon=143.4;
gridflag=1;
flagsource=1;
%}

%ADD PATHS FOR FUNCTIONS
addpath ../script_utils;
addpath ./gridgen;


%RUN THE GRID GENERATION CODE
if (~exist('xg'))
%    xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);   
%    xg=makegrid_tilteddipole_nonuniform_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);   
%    xg=makegrid_tilteddipole_nonuniform_oneside_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
  xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
%  xg=makegrid_cart_3D(xdist,lxp,ydist,lyp,I,glat,glon);
%  xg=makegrid_cart_3D_lowresx1(xdist,lxp,ydist,lyp,I,glat,glon);
end
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);


%IDENTIFICATION FOR THE NEW SIMULATION THAT IS TO BE DONE
simid='tohoku_medres'


%ALTERNATIVELY WE MAY WANT TO READ IN AN EXISTING OUTPUT FILE AND DO SOME INTERPOLATION ONTO A NEW GRID
fprintf('Reading in source file...\n');
%ID='~/simulations/3DPCarc_eq/'
%ID='~/zettergmdata/simulations/isinglass_eq/'
%ID='~/simulations/nepal2015_eq/'
%ID='~/zettergmdata/simulations/chile20153D_0.5_eq/'
%ID='~/zettergmdata/simulations/2Dtest_eq/'
ID='~/zettergmdata/simulations/tohoku_eq/'


%READ IN THE SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([ID,'/inputs/config.dat']);
xgin=readgrid(ID);
addpath ../vis/
direc=ID;


%FIND THE DATE OF THE END FRAEM OF THE SIMULATION (PRESUMABLY THIS WILL BE THE STARTING POITN FOR ANOTEHR)
[ymdend,UTsecend]=dateinc(tdur,ymd0,UTsec0);
filename=datelab(ymdend,UTsecend);
filename=[filename,'.dat']

%filename='20170303_09000.000000.dat';


%LOAD THE FRAME
loadframe3Dcurv;
rmpath ../vis/


%DO THE INTERPOLATION
if (lx3~=1)
  fprintf('Starting interp3''s...\n');
  [X2,X1,X3]=meshgrid(xgin.x2(3:end-2),xgin.x1(3:end-2),xgin.x3(3:end-2));
  [X2i,X1i,X3i]=meshgrid(xg.x2(3:end-2),xg.x1(3:end-2),xg.x3(3:end-2));
  for isp=1:lsp
    tmpvar=interp3(X2,X1,X3,ns(:,:,:,isp),X2i,X1i,X3i);
    inds=find(isnan(tmpvar));
    tmpvar(inds)=1e0;
    nsi(:,:,:,isp)=tmpvar;
    tmpvar=interp3(X2,X1,X3,vs1(:,:,:,isp),X2i,X1i,X3i);
    tmpvar(inds)=0e0;
    vs1i(:,:,:,isp)=tmpvar;
    tmpvar=interp3(X2,X1,X3,Ts(:,:,:,isp),X2i,X1i,X3i);
    tmpvar(inds)=100e0;
    Tsi(:,:,:,isp)=tmpvar;
  end
else
  fprintf('Starting interp2''s...\n');
  [X2,X1]=meshgrid(xgin.x2(3:end-2),xgin.x1(3:end-2));
  [X2i,X1i]=meshgrid(xg.x2(3:end-2),xg.x1(3:end-2));
  for isp=1:lsp
    tmpvar=interp2(X2,X1,squeeze(ns(:,:,:,isp)),X2i,X1i);
    inds=find(isnan(tmpvar));
    tmpvar(inds)=1e0;
    nsi(:,:,:,isp)=reshape(tmpvar,[lx1,lx2,1]);
    tmpvar=interp2(X2,X1,squeeze(vs1(:,:,:,isp)),X2i,X1i);
    tmpvar(inds)=0e0;
    vs1i(:,:,:,isp)=reshape(tmpvar,[lx1,lx2,1]);
    tmpvar=interp2(X2,X1,squeeze(Ts(:,:,:,isp)),X2i,X1i);
    tmpvar(inds)=100e0;
    Tsi(:,:,:,isp)=reshape(tmpvar,[lx1,lx2,1]);
  end
end


%WRITE OUT THE GRID
writegrid(xg,simid);    %just put it in pwd for now
dmy=[ymdend(3),ymdend(2),ymdend(1)];
writedata(dmy,UTsecend,nsi,vs1i,Tsi,simid);


%MAKE A SAMPLE PLOT OF INTERPOLATED DATA
%addpath ./vis;
%plotslice3D_curv(UTsec/3600,[3,11,2011],xgin,log10(ns(:,:,:,7)),'log_{10} n_e',[8 13])
%print -dpng ne.png
%plotslice3D_curv(UTsec/3600,[3,11,2011],xg,log10(nsi(:,:,:,7)),'log_{10} n_e',[8 13])
%print -dpng nei.png
%rmpath ./vis;


%RESET PATHS
rmpath ../script_utils;
rmpath ./gridgen;


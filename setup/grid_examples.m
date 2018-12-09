%% ADD PATHS TO THE GRID GENERATION SCRIPTS
addpath ./gridgen;
addpath ../script_utils;


%% MOORE, OK GRID (FULL)
dtheta=25;
dphi=35;
lp=125;
lq=250;
lphi=40;
altmin=80e3;
glat=39;
glon=262.51;
gridflag=0;
flagsource=1;


%% GEOGRAPHIC COORDINATES OF NEUTRAL SOURCE (OR GRID CENTER)
%MOORE OK
sourcelat=35.3;
sourcelong=360-97.7;
neugridtype=0;            %1 = Cartesian neutral grid, anything else - axisymmetric
zmin=0;
zmax=660;
rhomax=1800;

% % NO SOURCE SPECIFIED, SET TO CENTER OF GRID
% sourcelat=glat;
% sourcelong=glon;


%% FOR USERS INFO CONVERT SOURCE LOCATION TO GEOMAG
[sourcetheta,sourcephi]=geog2geomag(sourcelat,sourcelong);
sourcemlat=90-sourcetheta*180/pi;
sourcemlon=sourcephi*180/pi;


%% RUN THE GRID GENERATION CODE
if (~exist('xg'))
    xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);   
%    xg=makegrid_tilteddipole_nonuniform_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);   
%    xg=makegrid_tilteddipole_nonuniform_oneside_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
%     xg=makegrid_cart_3D(xdist,lxp,ydist,lyp,I,glat,glon);
end


ha=plotgrid(xg,flagsource,sourcelat,sourcelong,neugridtype,zmin,zmax,rhomax);


%% RETURN PATH VARIABLES TO NORMAL
rmpath ./gridgen;
rmpath ../script_utils;


%% ADDITIONAL EXAMPLES OF GRIDS AND SOURCE LOCATIONS...



% %TOHOKU
% sourcelat=38.429575;
% sourcelong=142.734757;

% %CHILE 2015
% sourcelat=-31.57;
% sourcelong=360-71.654;

% %NEPAL 2015
% sourcelat=28.17;
% sourcelong=85.48;

% %CHILE 2010
% sourcelat=-36.122;
% sourcelong=360-72.898;



% %TOHOKU-LIKE GRID
% dtheta=7.5;
% dphi=8.5;    %to make sure we encapsulate the neutral grid long. extent
% lp=50;
% lq=350;
% lphi=30;
% altmin=80e3;
% %glat=43.95;    %WRONG!!!
% glat=42.45;
% glon=143.4;
% gridflag=1;

% %GEOGRAPHICALLY CORRECT TOHOKU GRID
% dtheta=9;
% dphi=15;    %to make sure we encapsulate the neutral grid long. extent
% lp=50;
% %lp=350;    %4 km resolution
% %lq=150;
% lq=1000;
% %lq=1500;    %3-5km grid resolution
% lphi=30;
% %lphi=350;   %3-5km res.
% altmin=80e3;
% glat=42.45;
% glon=143.4;
% gridflag=1;

% %A THINNER TOHOKU GRID
% dtheta=7.5;
% dphi=12;
% lp=50;
% lq=250;
% lphi=50;
% altmin=80e3;
% glat=42.45;
% glon=143.4;
% gridflag=1;

% %CHILE 2015 GRID
% dtheta=8;
% dphi=14;
% lp=50;
% lq=250;
% lphi=50;
% altmin=80e3;
% glat=17.0;
% glon=288.2;
% gridflag=1;

% %NEPAL 2015 GRID
% dtheta=10;
% dphi=16;
% lp=100;
% lq=350;
% lphi=100;
% altmin=80e3;
% glat=35.75;
% glon=84.73;
% gridflag=1;

% %MOORE, OK GRID (PARTIAL)
% dtheta=15;
% dphi=20;
% lp=150;
% lq=500;
% lphi=25;
% altmin=80e3;
% glat=39;
% glon=262.51;
% gridflag=0;

% %RISR LOWRES GRID (CARTESIAN)
% xdist=1e6;
% ydist=1e6;
% lxp=50;
% lyp=50;
% glat=75.6975;
% glon=360.0-94.8322;
% gridflag=0;
% I=90;

% %CHILE 2010 GRID
% dtheta=8;
% dphi=14;
% lp=50;
% lq=250;
% lphi=50;
% altmin=80e3;
% glat=19.25;
% glon=287.5;
% gridflag=1;

% %CHILE 2010 GRID
% dtheta=8;
% dphi=14;
% lp=50;
% lq=250;
% lphi=50;
% altmin=80e3;
% glat=17.0;
% glon=288.2;
% gridflag=1;

% %PFISR LOWRES GRID (CARTESIAN)
% xdist=640e3;    %eastward distance
% ydist=385e3;    %northward distance
% lxp=100;
% lyp=100;
% glat=67.11;
% glon=212.95;
% gridflag=0;
% I=90;

% %EXAMPLE FOR KRISTINA MIDEX MISSION
% xdist=600e3;    %eastward distance
% ydist=100e3;    %northward distance
% lxp=384;
% lyp=512;
% glat=67.11;
% glon=212.95;
% gridflag=0;
% I=90;
% flagsource=0;     %specify no source

% %A HIGHRES TOHOKU
% dtheta=7.5;
% dphi=12;
% lp=50;
% lq=250;
% lphi=50;
% altmin=80e3;
% glat=42.45;
% glon=143.4;
% gridflag=1;
% flagsource=1;

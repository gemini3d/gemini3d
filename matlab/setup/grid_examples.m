function grid_examples()

%% Iowa grid for AGU 2019
dtheta=16;
dphi=29;
lp=100;
lq=200;
lphi=40;
altmin=80e3;
%glat=40;   %38.9609;
glat=41.5;   %38.9609;
glon=360-94.088;
gridflag=1;
flagsource=1;
iscurv=true;


%% GEOGRAPHIC COORDINATES OF NEUTRAL SOURCE (OR GRID CENTER)
% Iowa example
neuinfo.sourcelat=38.9609;
neuinfo.sourcelong=360-94.088;
neuinfo.neugridtype=3;    %1 = Cartesian neutral grid (2D), 2 - axisymmetric (2D), 3 - 3D Cartesian
neuinfo.zmin=0;
neuinfo.zmax=375;
neuinfo.xmin=-1200;
neuinfo.xmax=1200;
neuinfo.ymin=-1200;
neuinfo.ymax=1200;
neuinfo.rhomax=[];        %meaningless in 3D situations


%% FOR USERS INFO CONVERT SOURCE LOCATION TO GEOMAG
[sourcetheta,sourcephi]=geog2geomag(neuinfo.sourcelat,neuinfo.sourcelong);
sourcemlat=90-sourcetheta*180/pi;
sourcemlon=sourcephi*180/pi;


%% RUN THE GRID GENERATION CODE
if ~exist('xg', 'var')
  if iscurv
    xg = makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
%    xg=makegrid_tilteddipole_nonuniform_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
%    xg=makegrid_tilteddipole_nonuniform_oneside_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
  else
    xg = makegrid_cart_3D(p);
  end
end


%% PLOT THE GRID AND NEUTRAL INPUT EXTENT
%ha=plotgrid(xg,flagsource,sourcelat,sourcelong,neugridtype,zmin,zmax,rhomax);
ha = plot_mapgrid(xg,flagsource,neuinfo);
end % function

%% ADDITIONAL EXAMPLES OF GRIDS AND SOURCE LOCATIONS...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Iowa example
% neuinfo.sourcelat=38.9609;
% neuinfo.sourcelong=360-94.088;
% neuinfo.neugridtype=3;    %1 = Cartesian neutral grid (2D), 2 - axisymmetric (2D), 3 - 3D Cartesian
% neuinfo.zmin=0;
% neuinfo.zmax=375;
% neuinfo.xmin=-1200;
% neuinfo.xmax=1200;
% neuinfo.ymin=-1200;
% neuinfo.ymax=1200;
% neuinfo.rhomax=[];        %meaningless in 3D situations

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


% %MOORE OK
% sourcelat=35.3;
% sourcelong=360-97.7;
% neugridtype=0;            %1 = Cartesian neutral grid, anything else - axisymmetric
% zmin=0;
% zmax=660;
% rhomax=1800;

% % NO SOURCE
% sourcelat=glat;
% sourcelong=glon;
% neugridtype=[];               %1 = Cartesian neutral grid, anything else - axisymmetric
% zmin=[];
% zmax=[];
% rhomax=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% MSTIDs grid, CONUS
% dtheta=20;
% dphi=27.5;
% lp=128;
% lq=256;
% lphi=40;
% altmin=80e3;
% glat=39;
% glon=262.51;
% gridflag=0;
% flagsource=0;
% iscurv=true;

% %% MOORE, OK GRID (FULL)
% dtheta=25;
% dphi=35;
% lp=125;
% lq=250;
% lphi=40;
% altmin=80e3;
% glat=39;
% glon=262.51;
% gridflag=0;
% flagsource=1;

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

% %% EXAMPLE GRID CENTERED ON LONGYEAR
% dtheta=15;
% dphi=15;
% lp=125;
% lq=250;
% lphi=40;
% altmin=80e3;
% glat=78;
% glon=15;
% gridflag=0;
% flagsource=0;
% iscurv=true;

% %EXAMPLE FOR KRISTINA MIDEX MISSION
% xdist=600e3;    %eastward distance
% ydist=600e3;    %northward distance
% lxp=128;
% lyp=128;
% glat=67.11;
% glon=212.95;
% gridflag=0;
% I=90;
% flagsource=0;     %specify no source
% iscurv=false;

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
% flagsource=0;
% iscurv=true;

% %SAPs grid
% dtheta=15;
% dphi=75;
% lp=128;
% lq=256;
% lphi=64;
% altmin=80e3;
% glat=45;
% glon=262.51;
% gridflag=0;
% flagsource=0;
% iscurv=true;

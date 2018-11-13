function ha=plotgrid(xg,flagsource,sourcelat,sourcelong,neugridtype,zmin,zmax,rhomax)


%% INPUT COORDS NEED TO BE CONVERTED TO MAGNETIC
[sourcetheta,sourcephi]=geog2geomag(sourcelat,sourcelong);
sourcemlat=90-sourcetheta*180/pi;
sourcemlon=sourcephi*180/pi;
if (360-sourcemlon<20)
    sourcemlonplot=sourcemlon-360;
else
    sourcemlonplot=sourcemlon;
end


%% SET UP A MAP PLOT IF MAPPING TOOLBOX EXISTS
figure; hold on;
if (license('test','Map_Toolbox'))
  ha=gca;
  axesm('MapProjection','Mercator','MapLatLimit',[-(abs(sourcemlat)+30),abs(sourcemlat)+30],'MapLonLimit',[sourcemlonplot-40,sourcemlonplot+40])
  plotfun=@plot3m;
else
  plotfun=@plot3;   %no mapping toolbox so we'll do normal plots - the annoying thing is that these should be lon,lat,alt
end


%% CONVERT INPUT GRID COORDINATES INTO MLAT,MLON,ALT
Re=6370e3;
mlat=90-xg.theta*180/pi;
dphi=max(xg.phi(:))-min(xg.phi(:));
mlon=xg.phi*180/pi;
alt=xg.alt/1e3;
altscale=max(alt(:));
alt=alt/altscale;

if (360-sourcemlon<20)    %shift coords. too close to edge
   inds=find(mlon>180);
   mlon(inds)=mlon(inds)-360;
end


%% PLOT THE OUTLINE OF THE GRID
LW=2;
altlinestyle=':';    %line style for the "back" of the grid

h=plotfun(mlat(:,1,1),mlon(:,1,1),alt(:,1,1),'LineWidth',LW);
plotfun(mlat(:,1,end),mlon(:,1,end),alt(:,1,end),altlinestyle,'LineWidth',LW);
plotfun(mlat(:,end,1),mlon(:,end,1),alt(:,end,1),'LineWidth',LW);
h=plotfun(mlat(:,end,end),mlon(:,end,end),alt(:,end,end),altlinestyle,'LineWidth',LW);
linecolor=h.Color;

x=squeeze(mlat(1,:,1));
y=squeeze(mlon(1,:,1));
z=squeeze(alt(1,:,1));
plotfun(x,y,z,'LineWidth',LW);

x=squeeze(mlat(1,:,end));
y=squeeze(mlon(1,:,end));
z=squeeze(alt(1,:,end));
plotfun(x,y,z,altlinestyle,'LineWidth',LW);

x=squeeze(mlat(end,:,1));
y=squeeze(mlon(end,:,1));
z=squeeze(alt(end,:,1));
plotfun(x,y,z,'LineWidth',LW);

x=squeeze(mlat(end,:,end));
y=squeeze(mlon(end,:,end));
z=squeeze(alt(end,:,end));
plotfun(x,y,z,altlinestyle,'LineWidth',LW);

ix3=floor(xg.lx(3)/2);
x=squeeze(mlat(1,1,1:ix3));
y=squeeze(mlon(1,1,1:ix3));
z=squeeze(alt(1,1,1:ix3));
plotfun(x,y,z,'LineWidth',LW);

x=squeeze(mlat(1,1,ix3:xg.lx(3)));
y=squeeze(mlon(1,1,ix3:xg.lx(3)));
z=squeeze(alt(1,1,ix3:xg.lx(3)));
plotfun(x,y,z,altlinestyle,'LineWidth',LW);

x=squeeze(mlat(1,end,1:ix3));
y=squeeze(mlon(1,end,1:ix3));
z=squeeze(alt(1,end,1:ix3));
plotfun(x,y,z,'LineWidth',LW);

x=squeeze(mlat(1,end,ix3:xg.lx(3)));
y=squeeze(mlon(1,end,ix3:xg.lx(3)));
z=squeeze(alt(1,end,ix3:xg.lx(3)));
plotfun(x,y,z,altlinestyle,'LineWidth',LW);

x=squeeze(mlat(end,1,1:ix3));
y=squeeze(mlon(end,1,1:ix3));
z=squeeze(alt(end,1,1:ix3));
plotfun(x,y,z,'LineWidth',LW);

x=squeeze(mlat(end,1,ix3:xg.lx(3)));
y=squeeze(mlon(end,1,ix3:xg.lx(3)));
z=squeeze(alt(end,1,ix3:xg.lx(3)));
plotfun(x,y,z,altlinestyle,'LineWidth',LW);

x=squeeze(mlat(end,end,1:ix3));
y=squeeze(mlon(end,end,1:ix3));
z=squeeze(alt(end,end,1:ix3));
plotfun(x,y,z,'LineWidth',LW);

x=squeeze(mlat(end,end,ix3:xg.lx(3)));
y=squeeze(mlon(end,end,ix3:xg.lx(3)));
z=squeeze(alt(end,end,ix3:xg.lx(3)));
plotfun(x,y,z,altlinestyle,'LineWidth',LW);

xlabel('magnetic longitude (deg.)');
ylabel('magnetic latitude (deg.)');
zlabel('altitidue (km)');

h=gca;     %now go back and make all lines the same color...
lline=numel(h.Children);
for iline=1:16    %the last three children are the surface and text label objects
    h.Children(iline).Color=[0 0 0];
end


%% PLOT THE NEUTRAL SOURCE LOCATION, IF REQUESTED
if (flagsource~=0)
    hold on;
    plotfun(sourcemlat,sourcemlonplot,0,'ro','MarkerSize',16);
    
    
    %NOW CREATE A NEUTRAL GRID AND OVERPLOT IT
    %zmin=0;
    %%zmax=750;    %most earthquakes
    %zmax=660;    %Moore, OK
    lz=750;
    rhomin=0;
    %%rhomax=750;    %most earthquakes
    %rhomax=1800;    %Moore, OK
    lrho=750;
    zn=linspace(zmin,zmax,lz);
    rhon=linspace(rhomin,rhomax,lrho);
    rn=6370+zn;    %geocentric distance (in km)
    
    drho=rhomax-rhomin;                                                  %radius of circle, in kilometers, describing perp. directions of axisymmetric model
    xn=linspace(-1*drho,drho,100);                                       %N-S distance spanned by neutral model ("fake" number of grid points used here)
    dthetan=(max(xn(:))-min(xn(:)))/rn(1);                               %equivalent theta coordinates of the neutral mesh (used in the plot of grid)
    thetan=linspace(sourcetheta-dthetan/2,sourcetheta+dthetan/2,100);    %theta coordinates of N-S distance specified
    phinhalf1=sourcephi+sqrt((dthetan/2)^2-(thetan-sourcetheta).^2);
    phinhalf2=sourcephi-sqrt((dthetan/2)^2-(thetan-sourcetheta).^2);
    mlatn=90-thetan*180/pi;
    mlonnhalf1=phinhalf1*180/pi;
    mlonnhalf2=phinhalf2*180/pi;
    linemlon=sourcephi*ones(size(zn))*180/pi;
    
    
    %CORRECT THE MLONS TOO CLOSE TO EDGE OF GRID
    if (360-sourcemlon<20)
        inds=find(mlonnhalf1>180);
        mlonnhalf1(inds)=mlonnhalf1(inds)-360;
        inds=find(mlonnhalf2>180);
        mlonnhalf2(inds)=mlonnhalf2(inds)-360;
        linemlon=linemlon-360;
    end
    
    
    
    %PLOT THE NEUTRAL GRID
    hold on;
    zn=zn/altscale;
    plotfun(mlatn,real(mlonnhalf1),zn(1)*ones(size(mlatn)),altlinestyle,'LineWidth',LW);
    plotfun(mlatn,real(mlonnhalf2),zn(1)*ones(size(mlatn)),'LineWidth',LW);
    plotfun(mlatn,real(mlonnhalf1),zn(end)*ones(size(mlatn)),altlinestyle,'LineWidth',LW);
    plotfun(mlatn,real(mlonnhalf2),zn(end)*ones(size(mlatn)),'LineWidth',LW);
    plotfun(min(mlatn)*ones(size(zn)),linemlon,zn,'LineWidth',LW);
    plotfun(max(mlatn)*ones(size(zn)),linemlon,zn,'LineWidth',LW);
    hold off;
    
    
    %GO BACK AND MAKE ALL NEUTRAL GRID LINES THE SAME COLOR
    h=gca;
    lline=numel(h.Children);
    for iline=1:6    %the line objects for each axis are added in a "stack" fashion (last in first out)
        h.Children(iline).Color=linecolor;
    end
end


%CLEAN UP FONT SIZES, ETC.
FS=12;
set(gca,'FontSize',FS);
grid on;
set(gca,'LineWidth',LW-0.5)
view(1,2);


%CARTESIAN NEUTRAL GRID
%{
dz=72e3;
drho=72e3;

lz=9;
lrho=6;

zn=linspace(0,lz*dz,lz)';
rhon=linspace(0,lrho*drho,lrho);
xn=[-1*fliplr(rhon),rhon(2:lrho)];
lx=numel(xn);
yn=xn;          %all based off of axisymmetric model
rn=zn+6370e3;   %convert altitude to geocentric distance

dtheta=(max(xn(:))-min(xn(:)))/rn(1);    %equivalent theta coordinates of the neutral mesh
dphi=(max(yn(:))-min(yn(:)))/rn(1);      %should be a sin(theta) there?
thetan=linspace(meanth-dtheta/2,meanth+dtheta/2,2*lx-1);    %fudge this to make it look good.
phin=linspace(meanphi-dphi/2,meanphi+dphi/2,2*lx-1);
[THETAn,PHIn,Rn]=meshgrid(thetan,phin,rn);

MLATn=90-THETAn*180/pi;
MLONn=PHIn*180/pi;
Zn=(Rn-6370e3)/1e3;

hold on;
plot3(MLATn(:,end,1),MLONn(:,end,1),Zn(:,end,1),'LineWidth',LW);
h=plot3(MLATn(:,end,end),MLONn(:,end,end),Zn(:,end,end),'LineWidth',LW);
linecolor=h.Color;    %grab the color of the second line
plot3(MLATn(:,1,1),MLONn(:,1,1),Zn(:,1,1),'LineWidth',LW);
plot3(MLATn(:,1,end),MLONn(:,1,end),Zn(:,1,end),'LineWidth',LW);
plot3(squeeze(MLATn(1,:,1)),squeeze(MLONn(1,:,1)),squeeze(Zn(1,:,1)),'LineWidth',LW);
plot3(squeeze(MLATn(1,:,end)),squeeze(MLONn(1,:,end)),squeeze(Zn(1,:,end)),'LineWidth',LW);
plot3(squeeze(MLATn(end,:,1)),squeeze(MLONn(end,:,1)),squeeze(Zn(end,:,1)),altlinestyle,'LineWidth',LW);
plot3(squeeze(MLATn(end,:,end)),squeeze(MLONn(end,:,end)),squeeze(Zn(end,:,end)),altlinestyle,'LineWidth',LW);
plot3(squeeze(MLATn(1,1,:)),squeeze(MLONn(1,1,:)),squeeze(Zn(1,1,:)),'LineWidth',LW);
plot3(squeeze(MLATn(1,end,:)),squeeze(MLONn(1,end,:)),squeeze(Zn(1,end,:)),'LineWidth',LW);
plot3(squeeze(MLATn(end,1,:)),squeeze(MLONn(end,1,:)),squeeze(Zn(end,1,:)),altlinestyle,'LineWidth',LW);
plot3(squeeze(MLATn(end,end,:)),squeeze(MLONn(end,end,:)),squeeze(Zn(end,end,:)),altlinestyle,'LineWidth',LW);
%}


%% ADD MAPPED COASTLINES TO PLOT - DONE LAST TO MAKE SOME OF THE GRID LINE PLOTTING CODE CLEANER
if (license('test','Map_Toolbox'))
  hold on;
  ax=axis;
  load coastlines;
  [thetacoast,phicoast]=geog2geomag(coastlat,coastlon);
  mlatcoast=90-thetacoast*180/pi;
  mloncoast=phicoast*180/pi;
  
  if (360-sourcemlon<20)
    inds=find(mloncoast>180);
    mloncoast(inds)=mloncoast(inds)-360;
  end
  
  plotfun(mlatcoast,mloncoast,zeros(size(mlatcoast)),'b-','LineWidth',0.5);
  setm(gca,'MeridianLabel','on','ParallelLabel','on','MLineLocation',10,'PLineLocation',10,'MLabelLocation',10,'PLabelLocation',10);
  hold off;
end
view(270,35);
axis tight;


% %MAKE A MOVIE OF THE GRID ROTATING
% direc='~/Downloads/gridplot/';
% system(['mkdir ',direc])
% azstart=255;
% az=azstart:1:azstart+359;
% el=35;
% for iaz=1:numel(az)
%     view(az(iaz),el);
%     azstr=num2str(az(iaz));
%     ndigits=floor(log10(az(iaz)));
%     while(ndigits<2)
%        azstr=['0',azstr]; 
%        ndigits=ndigits+1;
%     end
%     print('-dpng',[direc,azstr,'.png'],'-r300');
%     %print('-depsc2',[direc,azstr,'.eps']);    
% end


end %function plotgrid
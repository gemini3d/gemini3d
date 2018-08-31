%EXAMPLE PARAMETERS
dtheta=10;
dphi=10;
lp=35;
lq=200;
lphi=48;
altmin=80e3;
glat=25;    %low-latitude
glon=270;
gridflag=1;

% %EXAMPLE PARAMETERS
% dtheta=10;
% dphi=10;
% lp=35;
% lq=350;
% lphi=48;
% altmin=80e3;
% glat=25;    %low-latitude
% glon=270;
% gridflag=1;

% dtheta=2;
% dphi=5;
% lp=35;
% lq=200;
% lphi=44;
% altmin=80e3;
% glat=65;    %high-latitude
% glon=270;
% gridflag=0;


%RUN THE GRID GENERATION CODE
xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);
Re=6370e3;


%TEST THE RESOLUTION IN THE DIFFERENT DIRECTIONS
ix3=floor(xg.lx(3)/2);    %can use a fixed x3 location since metric factors don't change with x3
ix2=1;
dl11=xg.h1(1:end-1,ix2,ix3).*xg.dx1b;
ix2=floor(xg.lx(2)/2);
dl12=xg.h1(1:end-1,ix2,ix3).*xg.dx1b;
ix2=xg.lx(2);
dl13=xg.h1(1:end-1,ix2,ix3).*xg.dx1b;

ix1=1;
dl21=squeeze(xg.h2(ix1,1:end-1,ix3)).*squeeze(xg.dx2b)';
ix1=floor(xg.lx(1)/2);
dl22=squeeze(xg.h2(ix1,1:end-1,ix3)).*squeeze(xg.dx2b)';
ix1=xg.lx(1);
dl23=squeeze(xg.h2(ix1,1:end-1,ix3)).*squeeze(xg.dx2b)';

ix2=floor(xg.lx(2)/2);
ix1=1;
dl31=squeeze(xg.h3(ix1,ix2,1:end-1)).*squeeze(xg.dx3b);
ix1=floor(xg.lx(1)/2);
dl32=squeeze(xg.h3(ix1,ix2,1:end-1)).*squeeze(xg.dx3b);
ix1=xg.lx(1);
dl33=squeeze(xg.h3(ix1,ix2,1:end-1)).*squeeze(xg.dx3b);

figure;
subplot(131);
plot(dl11/1e3);
hold on;
plot(dl12/1e3);
plot(dl13/1e3);

subplot(132);
plot(dl21/1e3);
hold on;
plot(dl22/1e3);
plot(dl23/1e3);

subplot(133);
plot(dl31/1e3);
hold on;
plot(dl32/1e3);
plot(dl33/1e3);


%WRITE THE GRID DATA TO TO A BINARY FILE
%writegrid(xg,'');
%save -v7.3 curvtest.mat xg;   %to compare file sizes - seems MATLAB does some aggressive compression
%xg=readgrid('curvtest');


%PLOT THE OUTLINES OF THE DIPOLE GRID - THIS VERSION USES X,Y,Z COORDS.
figure;
hold on;

h=plot3(xg.x(:,1,1)/Re,xg.y(:,1,1)/Re,xg.z(:,1,1)/Re);
linecolor=h.Color;
plot3(xg.x(:,1,end)/Re,xg.y(:,1,end)/Re,xg.z(:,1,end)/Re,'--');
plot3(xg.x(:,end,1)/Re,xg.y(:,end,1)/Re,xg.z(:,end,1)/Re);
plot3(xg.x(:,end,end)/Re,xg.y(:,end,end)/Re,xg.z(:,end,end)/Re,'--');

x=squeeze(xg.x(1,:,1))/Re;
y=squeeze(xg.y(1,:,1))/Re;
z=squeeze(xg.z(1,:,1))/Re;
plot3(x,y,z);

x=squeeze(xg.x(1,:,end))/Re;
y=squeeze(xg.y(1,:,end))/Re;
z=squeeze(xg.z(1,:,end))/Re;
plot3(x,y,z,'--');

x=squeeze(xg.x(end,:,1))/Re;
y=squeeze(xg.y(end,:,1))/Re;
z=squeeze(xg.z(end,:,1))/Re;
plot3(x,y,z);

x=squeeze(xg.x(end,:,end))/Re;
y=squeeze(xg.y(end,:,end))/Re;
z=squeeze(xg.z(end,:,end))/Re;
plot3(x,y,z,'--');


x=squeeze(xg.x(1,1,:))/Re;
y=squeeze(xg.y(1,1,:))/Re;
z=squeeze(xg.z(1,1,:))/Re;
plot3(x,y,z);

x=squeeze(xg.x(1,end,:))/Re;
y=squeeze(xg.y(1,end,:))/Re;
z=squeeze(xg.z(1,end,:))/Re;
plot3(x,y,z);

x=squeeze(xg.x(end,1,:))/Re;
y=squeeze(xg.y(end,1,:))/Re;
z=squeeze(xg.z(end,1,:))/Re;
plot3(x,y,z);

x=squeeze(xg.x(end,end,:))/Re;
y=squeeze(xg.y(end,end,:))/Re;
z=squeeze(xg.z(end,end,:))/Re;
plot3(x,y,z);

%axis equal;
view(-22,6)


%GO BACK AND MAKE ALL LINES THE SAME COLOR
h=gca;
lline=numel(h.Children);
for iline=1:lline
    h.Children(iline).Color=linecolor;
end


% %PLOT THE UNIT VECTORS TO CHECK DIRECTION
% figure;
% i1s=1:5:size(xg.x,1);
% i2s=1:5:size(xg.x,2);
% i3s=1:5:size(xg.x,3);
% quiver3(xg.x(i1s,i2s,i3s),xg.y(i1s,i2s,i3s),xg.z(i1s,i2s,i3s),xg.e1(i1s,i2s,i3s,1),xg.e1(i1s,i2s,i3s,2),xg.e1(i1s,i2s,i3s,3))
% 


%PLOT THE OUTLINES OF THE DIPOLE GRID - THIS VERSION USES LAT,LON.,ALT. COORDS.
figure;
hold on;

h=plot3(xg.glat(:,1,1),xg.glon(:,1,1),xg.alt(:,1,1));
linecolor=h.Color;
plot3(xg.glat(:,1,end),xg.glon(:,1,end),xg.alt(:,1,end),'--');
plot3(xg.glat(:,end,1),xg.glon(:,end,1),xg.alt(:,end,1));
plot3(xg.glat(:,end,end),xg.glon(:,end,end),xg.alt(:,end,end),'--');

x=squeeze(xg.glat(1,:,1));
y=squeeze(xg.glon(1,:,1));
z=squeeze(xg.alt(1,:,1));
plot3(x,y,z);

x=squeeze(xg.glat(1,:,end));
y=squeeze(xg.glon(1,:,end));
z=squeeze(xg.alt(1,:,end));
plot3(x,y,z,'--');

x=squeeze(xg.glat(end,:,1));
y=squeeze(xg.glon(end,:,1));
z=squeeze(xg.alt(end,:,1));
plot3(x,y,z);

x=squeeze(xg.glat(end,:,end));
y=squeeze(xg.glon(end,:,end));
z=squeeze(xg.alt(end,:,end));
plot3(x,y,z,'--');


x=squeeze(xg.glat(1,1,:));
y=squeeze(xg.glon(1,1,:));
z=squeeze(xg.alt(1,1,:));
plot3(x,y,z);

x=squeeze(xg.glat(1,end,:));
y=squeeze(xg.glon(1,end,:));
z=squeeze(xg.alt(1,end,:));
plot3(x,y,z);

x=squeeze(xg.glat(end,1,:));
y=squeeze(xg.glon(end,1,:));
z=squeeze(xg.alt(end,1,:));
plot3(x,y,z);

x=squeeze(xg.glat(end,end,:));
y=squeeze(xg.glon(end,end,:));
z=squeeze(xg.alt(end,end,:));
plot3(x,y,z);

%axis equal;
view(-22,6)

xlabel('glat');
ylabel('glon');
zlabel('alt');

%GO BACK AND MAKE ALL LINES THE SAME COLOR
h=gca;
lline=numel(h.Children);
for iline=1:lline
    h.Children(iline).Color=linecolor;
end


%PLOT THE OUTLINES OF THE DIPOLE GRID - THIS VERSION USES MLAT,MLON.,ALT. COORDS.
mlat=90-xg.theta*180/pi;
mlon=xg.phi*180/pi;
alt=xg.alt/1e3;
figure;
hold on;

altlinestyle=':'

h=plot3(mlat(:,1,1),mlon(:,1,1),alt(:,1,1));
linecolor=h.Color;
plot3(mlat(:,1,end),mlon(:,1,end),alt(:,1,end),altlinestyle);
plot3(mlat(:,end,1),mlon(:,end,1),alt(:,end,1));
plot3(mlat(:,end,end),mlon(:,end,end),alt(:,end,end),altlinestyle);

x=squeeze(mlat(1,:,1));
y=squeeze(mlon(1,:,1));
z=squeeze(alt(1,:,1));
plot3(x,y,z);

x=squeeze(mlat(1,:,end));
y=squeeze(mlon(1,:,end));
z=squeeze(alt(1,:,end));
plot3(x,y,z,altlinestyle);

x=squeeze(mlat(end,:,1));
y=squeeze(mlon(end,:,1));
z=squeeze(alt(end,:,1));
plot3(x,y,z);

x=squeeze(mlat(end,:,end));
y=squeeze(mlon(end,:,end));
z=squeeze(alt(end,:,end));
plot3(x,y,z,altlinestyle);


x=squeeze(mlat(1,1,:));
y=squeeze(mlon(1,1,:));
z=squeeze(alt(1,1,:));
plot3(x,y,z);

x=squeeze(mlat(1,end,:));
y=squeeze(mlon(1,end,:));
z=squeeze(alt(1,end,:));
plot3(x,y,z);

x=squeeze(mlat(end,1,:));
y=squeeze(mlon(end,1,:));
z=squeeze(alt(end,1,:));
plot3(x,y,z);

x=squeeze(mlat(end,end,:));
y=squeeze(mlon(end,end,:));
z=squeeze(alt(end,end,:));
plot3(x,y,z);

%axis equal;
view(-22,6)

xlabel('magnetic latitude (deg.)');
ylabel('magnetic longitude (deg.)');
zlabel('altitidue (km)');

%GO BACK AND MAKE ALL LINES THE SAME COLOR
h=gca;
lline=numel(h.Children);
for iline=1:lline
    h.Children(iline).Color=linecolor;
end


%NOW CREATE A NEUTRAL GRID AND OVERPLOT IT
meanth=mean(xg.theta(xg.lx(1),:,floor(end/2)))+5.75*pi/180;   %N lowlat, TD model
meanphi=mean(xg.phi(xg.lx(1),floor(end/2),:));
dcoord=74e3;    %resolution of the neutral dynamics model

lz=9;
lrho=6;

zn=linspace(0,lz*dcoord,lz)';
rhon=linspace(0,lrho*dcoord,lrho);
xn=[-1*fliplr(rhon),rhon(2:lrho)];
lx=numel(xn);
yn=xn;    %all based off of axisymmetric model
rn=zn+6370e3;   %convert altitude to geocentric distance

dtheta=(max(xn(:))-min(xn(:)))/rn(1);    %equivalent theta coordinates of the neutral mesh
dphi=(max(yn(:))-min(yn(:)))/rn(1);    %should be a sin(theta) there?
thetan=linspace(meanth-dtheta/2,meanth+dtheta/2,2*lx-1);
phin=linspace(meanphi-dphi/2,meanphi+dphi/2,2*lx-1);
[THETAn,PHIn,Rn]=meshgrid(thetan,phin,rn);

MLATn=90-THETAn*180/pi;
MLONn=PHIn*180/pi;
Zn=(Rn-6370e3)/1e3;

hold on;
plot3(MLATn(:,end,1),MLONn(:,end,1),Zn(:,end,1));
h=plot3(MLATn(:,end,end),MLONn(:,end,end),Zn(:,end,end));
linecolor=h.Color;    %grab the color of the second line
plot3(MLATn(:,1,1),MLONn(:,1,1),Zn(:,1,1));
plot3(MLATn(:,1,end),MLONn(:,1,end),Zn(:,1,end));
plot3(squeeze(MLATn(1,:,1)),squeeze(MLONn(1,:,1)),squeeze(Zn(1,:,1)));
plot3(squeeze(MLATn(1,:,end)),squeeze(MLONn(1,:,end)),squeeze(Zn(1,:,end)));
plot3(squeeze(MLATn(end,:,1)),squeeze(MLONn(end,:,1)),squeeze(Zn(end,:,1)),altlinestyle);
plot3(squeeze(MLATn(end,:,end)),squeeze(MLONn(end,:,end)),squeeze(Zn(end,:,end)),altlinestyle);
plot3(squeeze(MLATn(1,1,:)),squeeze(MLONn(1,1,:)),squeeze(Zn(1,1,:)));
plot3(squeeze(MLATn(1,end,:)),squeeze(MLONn(1,end,:)),squeeze(Zn(1,end,:)));
plot3(squeeze(MLATn(end,1,:)),squeeze(MLONn(end,1,:)),squeeze(Zn(end,1,:)),altlinestyle);
plot3(squeeze(MLATn(end,end,:)),squeeze(MLONn(end,end,:)),squeeze(Zn(end,end,:)),altlinestyle);
hold off;

%GO BACK AND MAKE ALL NEUTRAL GRID LINES THE SAME COLOR
h=gca;
lline=numel(h.Children);
for iline=1:lline/2    %the line objects for each axis are added in a "stack" fashion (last in first out)
    h.Children(iline).Color=linecolor;
%    h.Children(iline).LineWidth=0.60;
end

FS=18;
set(gca,'FontSize',FS);
view(-11,4);
grid on;
print -depsc2 3Datmosionos_grid.eps;

close all;
clear;
clc;
direc='~/gemini.fort.ind/output/'

filename='0.000001.dat'
loadframe3D;


%BEGIN YISHI'S PLOT CODE
ny = numel(x3);nx = numel(x2);nz = numel(x1);
exam_alt = [120;300;750];
exam_y = 0; 
exam_x = 0;

%Create Quiver Plots
num_arrows_x = 8;
num_arrows_y = 8;
num_arrows_z = 6;
locx = 10:floor((nx-20)/num_arrows_x):nx-10;
locy = 10:floor((ny-20)/num_arrows_y):ny-10;
locz = 10:floor((nz-20)/num_arrows_z):nz-10;
[yy,xx] = meshgrid(x3(locy),x2(locx));
[yy_z,zz] = meshgrid(x3(locy),x1(locz));
[xx_z,zzz] = meshgrid(x2(locx),x1(locz));


[~,alt_index1] = min(abs(x1-exam_alt(1)*1e3));
[~,alt_index2] = min(abs(x1-exam_alt(2)*1e3));
[~,alt_index3] = min(abs(x1-exam_alt(3)*1e3));
[~,y_index1] = min(abs(x3-exam_y(1)*1e3));
[~,x_index1] = min(abs(x2-exam_x(1)*1e3));


%SELECT THE PARAMETER FOR PLOTTING
param=ne;
cbar_title='\bf n_{e} [m^{-3}]';
%param=log10(ne);
%cbar_title='\bf log_{10} n_{e} [m^{-3}]';
%param=1e6*Jpar;
%cbar_title='\bf J_{||} [uA m^{-2}]';
%param=1e6*Jperp2;
%cbar_title='\bf J_{\perp,2} [uA m^{-2}]';
%param=1e6*Jperp3;
%cbar_title='\bf J_{\perp,3} [uA m^{-2}]';
%param=Ti;
%cbar_title='\bf T_{i} [K]';
%param=p;
%cbar_title='\bf n_{O^+}/n_e';
%param=vi;
%cbar_title='\bf v_{||} [m/s]';
%param=Te;
%cbar_title='\bf T_{e} [K]';


%MAKE THE PLOT
%fullscreen = get(0,'ScreenSize');
%h1=figure ('PaperType','A0','Position',[0 -50 fullscreen(3) fullscreen(4)]);
%set (h1, 'Units', 'normalized', 'Position', [0,0,1,1]);set(gcf,'PaperPosition',[0 0 14 8.5]);
figure;
set(gcf,'PaperPosition',[0 0 12 8.5]);


subplot(2,3,[2]); hold on;
s2 = [min(x3)/1e3;x3(y_index1)/1e3]; s1 = [min(x2)/1e3;x2(x_index1)/1e3]; s3 = [min(x1)/1e3]; %y,x,z
hsne = slice(X2/1e3,X3/1e3,X1/1e3,param,s1,s2,s3);
set(hsne,'FaceColor','interp','EdgeColor','none');axis tight
xlabel('\bf x [km]');
ylabel('\bf y [km]');
zlabel('\bf z [km]');
shading flat;
cax=caxis;
% cb=colorbar('westoutside');
% pos = get(cb,'pos');
% set(cb,'position',[0.05 pos(2) pos(3) pos(4)+0.44]);
% xlabel(cb,'\bf n_{e} [m^{-3}]');
ht=title(['\bf Simulation Time = ',num2str(t,'%.2f'),' sec'],'fontsize',14);
pos_tit = get(ht,'pos');
%set(ht,'position',[pos_tit(1) pos_tit(2) pos_tit(3)+600]);
set(ht,'position',[pos_tit(1) pos_tit(2) pos_tit(3)+1800]);
%view(95,3);
view(141,35)
axis equal;
hold off;

subplot(2,3,5); hold on;
s1 = min(x2)/1e3; s2 = min(x3)/1e3; s3 = [min(x1)/1e3;exam_alt]; %x,y,z
hsne = slice(X2/1e3,X3/1e3,X1/1e3,param,s1,s2,s3);
set(hsne,'FaceColor','interp','EdgeColor','none');axis tight
xlabel('\bf x [km]');
ylabel('\bf y [km]');
zlabel('\bf z [km]');
shading flat;
caxis(cax);
%cb=colorbar('westoutside');
%pos = get(cb,'pos');
%set(cb,'position',[0.37 pos(2) pos(3)-0.0083 pos(4)+0.44]);
%xlabel(cb,'\bf n_{e} [m^{-3}]');
%view(95,3);
view(141,25)
axis equal;
hold off;

subplot(2,3,[4]); 
pcolor(x2/1e3,x1/1e3,squeeze(param(y_index1,:,:))');
shading flat;%hc = colorbar;ylabel(hc,'\bf n_{e} [m^{-3}]'); 
caxis(cax);
cb=colorbar('westoutside');
pos = get(cb,'pos');
set(cb,'position',[0.07 pos(2) pos(3)-0.007 pos(4)+0.44]);
xlabel(cb,cbar_title);
hold on;
%h = quiver(xx_z/1e3,zzz/1e3,squeeze(vi2(y_index1,locx,locz))',squeeze(vi(y_index1,locx,locz))','color','w','autoscalefactor',1.3,'linewidth',1.5);
xlabel('\bf x [km]');
ylabel('\bf z [km]');
h = title(['\bf y = ',num2str(exam_y,'%d'),' km'],'position',[250 1010]);
hold off;

subplot(2,3,[1]); 
pcolor(x3/1e3,x1/1e3,squeeze(param(:,x_index1,:))');
shading flat;%hc = colorbar;ylabel(hc,'\bf n_{e} [m^{-3}]'); 
hold on;
caxis(cax);
%h = quiver(yy_z/1e3,zz/1e3,squeeze(vi3(locy,x_index1,locz))',squeeze(vi(locy,x_index1,locz))','color','w','autoscalefactor',1,'linewidth',1.5);
xlabel('\bf y [km]');
ylabel('\bf z [km]');
h = title(['\bf x = ',num2str(exam_x,'%d'),' km'],'position',[250 1010]);
hold off;

subplot(3,3,3);
pcolor(x2/1e3,x3/1e3,param(:,:,alt_index3));
shading flat;
cb=colorbar;
pos = get(cb,'pos');
set(cb,'position',[pos(1)+0.06 pos(2) pos(3)-0.007 pos(4)]);
ylabel(cb,cbar_title); 
hold on;
%h = quiver(xx/1e3,yy/1e3,squeeze(vi2(locy,locx,alt_index3))',squeeze(vi3(locy,locx,alt_index3))','color','w','autoscalefactor',1.1,'linewidth',1.5);
h = quiver(xx/1e3,yy/1e3,squeeze(vi2(locy,locx,alt_index3))',squeeze(vi3(locy,locx,alt_index3))','color','w','autoscalefactor',1.1,'linewidth',1.5);
xlabel('\bf x [km]');
ylabel('\bf y [km]');
%h = title(['\bf z = ',num2str(exam_alt(3),'%d'),' km'],'position',[250 1300]);
h = title(['\bf z = ',num2str(exam_alt(3),'%d'),' km'],'position',[250 1050]);
hold off;

subplot(3,3,6);
pcolor(x2/1e3,x3/1e3,param(:,:,alt_index2));
shading flat;
cb=colorbar;
pos = get(cb,'pos');
set(cb,'position',[pos(1)+0.06 pos(2) pos(3)-0.007 pos(4)]);
ylabel(cb,cbar_title); 
hold on;
%h = quiver(xx/1e3,yy/1e3,squeeze(vi2(locy,locx,alt_index2))',squeeze(vi3(locy,locx,alt_index2))','color','w','autoscalefactor',1.1,'linewidth',1.8);
h = quiver(xx/1e3,yy/1e3,squeeze(vi2(locy,locx,alt_index2))',squeeze(vi3(locy,locx,alt_index2))','color','w','autoscalefactor',1.1,'linewidth',1.8);
xlabel('\bf x [km]');
ylabel('\bf y [km]');
%title(['\bf z = ',num2str(exam_alt(2),'%d'),' km'],'position',[250 1300]);
h = title(['\bf z = ',num2str(exam_alt(2),'%d'),' km'],'position',[250 1050]);
hold off;

subplot(3,3,9);
pcolor(x2/1e3,x3/1e3,squeeze(param(:,:,alt_index1)));
shading flat;
cb=colorbar;
pos = get(cb,'pos');
set(cb,'position',[pos(1)+0.06 pos(2) pos(3)-0.007 pos(4)]);
ylabel(cb,cbar_title); 
hold on;
h = quiver(xx/1e3,yy/1e3,squeeze(vi2(locy,locx,alt_index1))',squeeze(vi3(locy,locx,alt_index1))','color','w','autoscalefactor',1.2,'linewidth',1.7);
xlabel('\bf x [km]');
ylabel('\bf y [km]');
%title(['\bf z = ',num2str(exam_alt(1),'%d'),' km'],'position',[800 470]);
h = title(['\bf z = ',num2str(exam_alt(1),'%d'),' km'],'position',[250 1050]);
hold off;









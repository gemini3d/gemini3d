close all;
clear;
clc;


%LOAD DATA FROM TIME OF INTEREST
%direc='~/simulations/2DPCarccommit/'
%direc='~/simulations/2DPCarcGLOW/'
%direc='~/simulations/GLOWGEMINI_eheating/'
%direc='~/simulations/2DGLOWtest/'
direc='/media/data/zettergm/simulations/gemini2D/2Ddiscarc/'
filename='20130220_18300.000000.dat'
loadframe3D;


%INTERPOLATION OF SIMULATION DATA ONTO A UNIFORM GRID
lsp=size(ns,4);
Re=6370e3;

izs=1:lx1;
ixs=1:lx2;
x=x3/1e3;
z=x1'/1e3;
lxp=500;     %number of points for interpolated grid in the x-direction
lzp=500;     %" for the z-direction
%minx=-125;
%maxx=125;
minx=-15;
maxx=15;
minz=90;     %min altitude to use in plots
maxz=500;    %max altitude to use in plots
xp=linspace(minx,maxx,lxp);    %use a uniformly spaced grid for interpolation so we can use imagesc instead of pcolor
zp=linspace(minz,maxz,lzp)';


%INTERPOLATE ONTO PLOTTING GRID
J1p=interp2(x,z,Jpar,xp,zp);
J2p=interp2(x,z,Jperp2,xp,zp);
nep=interp2(x,z,ne,xp,zp);
pp=interp2(x,z,p,xp,zp);
Tip=interp2(x,z,Ti,xp,zp);
Tep=interp2(x,z,Te,xp,zp);
vip=interp2(x,z,vi,xp,zp);
vi3p=interp2(x,z,vi3,xp,zp);


%COLORBAR LIMITES FOR EACH PARAMETER
vlims=[-max(abs(vip(:))) max(abs(vip(:)))];
v3lims=[-max(abs(vi3p(:))) max(abs(vi3p(:)))];
Tlims=[0 max(max(Tep(:)),max(Tip(:)))];
nlims=[0 max(nep(:))];    
Jlims=[-max(abs(J1p(:))/1e-6) max(abs(J1p(:))/1e-6)];
plab='n_{O^+}/n_e';
plims=[0 1];
FS=10;


%MAKE THE PLOTS!
figure;
set(gcf,'PaperPosition',[0 0 8.5 6]);


%CURRENTS IN FIRST ROW OF PLOTS

Jlab='J_1 (FAC) \n (u A m^{-2})';    %label for the colorbar of the top-left plot
Jlab2='J_3 (perp. B) \n (u A m^{-2})';    %label for the colorbar of the top-right plot
subplot(3,7,1:3);
h=imagesc(xp,zp,J1p*1e6);
set(gca,'FontSize',FS);
axis xy;
axis square;
axis tight;
caxis(Jlims)
c=colorbar;
xlabel(c,sprintf(Jlab))
xlabel('horizontal distance (km)');
ylabel('altitude (km)');
strval=datestr(datenum(simdate));
text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',10,'Color',[1 1 1]);

subplot(3,7,5:7);
h=imagesc(xp,zp,J2p*1e6);
set(gca,'FontSize',FS);
axis xy;
axis square;
axis tight;
%caxis(Jlims)
caxis([-0.5 0.5]);
c=colorbar;
xlabel(c,sprintf(Jlab2))
xlabel('horizontal distance (km)');
ylabel('altitude (km)');


%ION VELOCITY IN FIRST SET OF PLOTS
%{
vlab='v_1 (FAC) \n (m/s)';    %label for the colorbar of the top-left plot
vlab2='v_3 (perp. B) \n (m/s)';    %label for the colorbar of the top-right plot
subplot(3,7,1:3);
h=imagesc(xp,zp,vip);
set(gca,'FontSize',FS);
axis xy;
axis square;
axis tight;
caxis(vlims)
c=colorbar;
xlabel(c,sprintf(vlab))
xlabel('horizontal distance (km)');
ylabel('altitude (km)');
strval=datestr(datenum(simdate));
text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',10,'Color',[1 1 1]);

subplot(3,7,5:7);
h=imagesc(xp,zp,vi3p);
set(gca,'FontSize',FS);
axis xy;
axis square;
axis tight;
caxis(v3lims)
c=colorbar;
xlabel(c,sprintf(vlab2))
xlabel('horizontal distance (km)');
ylabel('altitude (km)');
%}

%REST OF PARAMETERS
subplot(3,7,8:10);
h=imagesc(xp,zp,nep);
set(gca,'FontSize',FS);
axis xy;
axis square;
axis tight;
caxis(nlims);
c=colorbar;
xlabel(c,sprintf('n_e (m^{-3})'))
xlabel('horizontal distance (km)');
ylabel('altitude (km)');

subplot(3,7,12:14);
h=imagesc(xp,zp,pp);
set(gca,'FontSize',FS);
axis xy;
axis square;
axis tight;
caxis(plims);
c=colorbar;
xlabel(c,plab)
xlabel('horizontal distance (km)');
ylabel('altitude (km)');

subplot(3,7,15:17);
h=imagesc(xp,zp,Tip);
set(gca,'FontSize',FS);
axis xy;
axis square;
axis tight;
caxis(Tlims)
c=colorbar;
xlabel(c,'T_i (K)')
xlabel('horizontal distance (km)');
ylabel('altitude (km)');

subplot(3,7,19:21);
h=imagesc(xp,zp,Tep);
set(gca,'FontSize',FS);
axis xy;
axis square;
axis tight;
caxis(Tlims)
c=colorbar;
xlabel(c,'T_e (K)')
xlabel('horizontal distance (km)');
ylabel('altitude (km)');

print -dpng -r300 frame.png;


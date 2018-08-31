close all;
clear;
clc;
%direc='/media/data/zettergm/simulations/preconditioned/'
direc='/media/data/zettergm/simulations/3DPCarc/'

filename='120.000000.dat'
loadframe3D;


%2D
%figure;
%subplot(131), pcolor(x3/1e3,x1/1e3,ne), shading flat, colormap(jet(256)), colorbar
%subplot(132), pcolor(x3/1e3,x1/1e3,Jpar), shading flat, colormap(jet(256)), colorbar
%subplot(133), pcolor(x3/1e3,x1/1e3,Ti), shading flat, colormap(jet(256)), colorbar

%3D
figure;
ix2=floor(lxs(2)/2);
subplot(131), pcolor(x3/1e3,x1/1e3,squeeze(ns(:,ix2,:,7))), shading flat, colormap(jet(256)), colorbar
subplot(132), pcolor(x3/1e3,x1/1e3,squeeze(J1(:,ix2,:))), shading flat, colormap(jet(256)), colorbar
subplot(133), pcolor(x3/1e3,x1/1e3,squeeze(Ts(:,ix2,:,1))), shading flat, colormap(jet(256)), colorbar


print -dpng test.png


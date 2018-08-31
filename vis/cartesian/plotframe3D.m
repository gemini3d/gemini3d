%FILE INFO
%direc='~/simulations/curvtest_tohoku_verylowres3/'
%filename='20110311_21863.000000.dat'
%direc='~/simulations/curvtest_tohoku_eq2/'
%filename='20110310_20783.000001.dat'
direc='~/simulations/curvtest_tohoku_highres_weak/'
filename='20110311_21443.000000.dat'


if (~exist('xg','var'))
  %LOAD THE DATA FROM THIS TIME
  %loadframe3Dcurv;
  loadframe3Dcurvavg;
  t=simdate(4);  %assumes fortran code outputs the date in UT hours.  

  %WE ALSO NEED TO LOAD THE GRID FILE
  xg=readgrid([direc,'/']);
  fprintf('Grid loaded...\n');
end


%ATTEMPT TO PLOT SOME INFO IN A LINE PLOT
%figure;
%h=plotprofiles3D(t/3600,xg,ns,vs1,Ts);
%print -dpng profiles.png;


%NOW ATTEMPT A 3D SLICE PLOT
% figure;
% h2=plotslice3D(t/3600,xg,Ts(:,:,:,7));
% figure;
% h2=plotslice3D(t/3600,xg,vs1(:,:,:,1));
% figure;
% h2=plotslice3D(t/3600,xg,log10(ns(:,:,:,7)));
% caxis([8 12])


%figure;
%h2=plotslice3D_curv_north(t/3600,[11,3,2011],xg,log10(ns(:,:,:,7)),'log n_e',[8 13]);
%[az,el]=view;
%view(az,10);
%print -dpng n.png;

%figure;
%h2=plotslice3D_curv_north(t/3600,[11,3,2011],xg,log10(ne),'log n_e',[8 12.7]);
%[az,el]=view;
%view(az,10);
%print -dpng n.png;

%figure;
%h2=plotslice3D_curv_north(t/3600,[11,3,2011],xg,Ts(:,:,:,7),'T_e',[0 5000]);
%[az,el]=view;
%view(az,10);
%print -dpng Te.png

%figure;
%h2=plotslice3D_curv_north(t/3600,[11,3,2011],xg,Te(:,:,:),'T_e',[0 5000]);
%[az,el]=view;
%view(az,10);
%print -dpng Te.png

%figure;
%h2=plotslice3D_curv_north(t/3600,[11,3,2011],xg,vs1(:,:,:,1),'v_1',[-350 350]);
%[az,el]=view;
%view(az,10);
%print -dpng v1.png

%{
figure;
h2=plotslice3D_curv_north(t,[11,3,2011],xg,v1(:,:,:),'v_1',[-100 100]);
[az,el]=view;
view(az,10);
print -dpng v1.png
%}

%{
figure;
h2=plot2D_curv_north(t,[11,3,2011],xg,v1(:,:,:),'v_1',[-300 300]);
print -dpng v1_2D.png
%}

%COMPUTE SOURUCE LOCATIOKN IN MCOORDS
glat=38.429575d0
glon=142.734757d0
addpath ../setup;
[theta,phi]=geog2geomag(glat,glon);
mlat=90-theta*180/pi;
mlon=phi*180/pi;
rmpath ../setup;

figure;
h2=plot2D_curv_north_frames(t,[11,3,2011],xg,v1(:,:,:),'v_1 (m/s)',[-100 100],[mlat,mlon]);
print -depsc2 v1_2D_frames.eps
print -dpng -r300 v1_2D_frames.png

figure;
h2=plot2D_curv_north_frames(t,[11,3,2011],xg,J1(:,:,:)*1e6,'J_1 (uA/m^2)',[-0.025 0.025],[mlat,mlon]);
print -depsc2 J1_2D_frames.eps
print -dpng -r300 J1_2D_frames.png

figure;
h2=plot2D_curv_north_frames(t,[11,3,2011],xg,v2(:,:,:),'v_2 (m/s)',[-5 5],[mlat,mlon]);
print -depsc2 v2_2D_frames.eps
print -dpng -r300 v2_2D_frames.png

figure;
h2=plot2D_curv_north_frames(t,[11,3,2011],xg,v3(:,:,:),'v_3 (m/s)',[-5 5],[mlat,mlon]);
print -depsc2 v3_2D_frames.eps
print -dpng -r300 v3_2D_frames.png

%{
figure;
h2=plotslice3D_curv_north_scattered(t,[11,3,2011],xg,v1(:,:,:),'v_1',[-100 100]);
[az,el]=view;
view(az,10);
print -dpng v1_scattered.png
%}


%{
figure;
h2=plotslice3D_curv_north_altinterp2(t,[11,3,2011],xg,v1(:,:,:),'v_1',[-100 100]);
[az,el]=view;
view(az,10);
print -dpng v1_slices.png
%}


%figure;
%h2=plotslice3D_curv_north(t/3600,[11,3,2011],xg,v2(:,:,:),'v_2',[-50 50]);
%[az,el]=view;
%view(az,10);
%print -dpng v2.png
%
%figure;
%h2=plotslice3D_curv_north(t/3600,[11,3,2011],xg,v3(:,:,:),'v_3',[-50 50]);
%[az,el]=view;
%view(az,10);
%print -dpng v3.png
%
%{
figure;
h2=plotslice3D_curv_north(t/3600,[11,3,2011],xg,J1(:,:,:)*1e6,'J_1',[-0.05 0.05]);
[az,el]=view;
view(az,10);
print -dpng J1.png

figure;
h2=plotslice3D_curv(t/3600,[11,3,2011],xg,J1(:,:,:)*1e6,'J_1',[-0.25 0.25]);
[az,el]=view;
view(az,10);
print -dpng J1full.png

figure;
minJ=min(J2(:))*1e6;
maxJ=max(J2(:))*1e6;
h2=plotslice3D_curv(t/3600,[11,3,2011],xg,J2(:,:,:)*1e6,'J_2',[minJ maxJ]);
[az,el]=view;
view(az,10);
print -dpng J2full.png

figure;
minJ=min(J3(:))*1e6;
maxJ=max(J3(:))*1e6;
h2=plotslice3D_curv(t/3600,[11,3,2011],xg,J3(:,:,:)*1e6,'J_3',[minJ maxJ]);
[az,el]=view;
view(az,10);
print -dpng J3full.png
%}

%{
figure;
h2=plotslice3D_curv(t/3600,xg,vs1(:,:,:,2));

figure;
h2=plotslice3D_curv(t/3600,xg,J1);
%}
% figure;
% h2=plotslice3D_curv_corner(t/3600,xg,log10(ns(:,:,:,7)));
% caxis([8 12])

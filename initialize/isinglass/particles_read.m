%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS PROGRAM REQUIRES THE RESTORE_IDL SCRIPTS WHICH CAN BE DOWNLOADED FROM THE
% MATLAB FILE EXCHANGE AT:  
%   https://www.mathworks.com/matlabcentral/fileexchange/43899-restore-idl 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idlpath='~/ISINGLASS/AGU2017/'
addpath([idlpath,'./restore_idl']);
addpath ../../script_utils;


%CREATE SOME SPACE FOR OUTPUT FILES
outdir='~/zettergmdata/simulations/input/isinglass_particles_highres/';
system(['mkdir ',outdir]);
system(['rm ',outdir,'/*']);   %clean out existing files


%READ IN THE IDL SAVE FILE - THIS IS THE FORMAT NORMALLY GIVEN TO ME BY GUY GRUBBS
%fname='isinglass_eflux_asi.sav';
%fname='isinglass_eflux_asi_sync.sav';
%fname='isinglass_eflux_asi_resync.sav';
datapath='~/ISINGLASS/AGU2017/';
fname='isinglass_eflux_asi_highres.sav';
outargs=restore_idl([datapath,fname]);
time=outargs.NEW_TIME;
lat=outargs.NEW_LAT;
lon=outargs.NEW_LON;
Qdat=outargs.RESAMP_Q;
E0dat=outargs.RESAMP_EO;
[lt,llon0,llat0]=size(Qdat);


%GEOGRAPHIC AND MAGNETIC POSITIONS OF DATA
glat=lat;
glon=lon;
gloncorrected=glon+360;    %for interpolations which won't understand periodic coordinate
% [theta,phi]=geog2geomag(glat,glon);
% mlat=90-theta*180/pi;
% mlon=phi*180/pi;


%CREATE A PLAID MLON,MLAT GRID FOR THE MODEL
theta=zeros(llon0,llat0);
phi=zeros(llon0,llat0);
for ilat=1:llat0
    for ilon=1:llon0
        [thetatmp,phitmp]=geog2geomag(glat(ilat),glon(ilon));
        theta(ilon,ilat)=thetatmp;
        phi(ilon,ilat)=phitmp;
    end
end

llat=250;
llon=250;
mlatdat=90-theta*180/pi;
mlondat=phi*180/pi;
mlatmin=min(mlatdat(:));
mlatmax=max(mlatdat(:));
mlat=linspace(mlatmin,mlatmax,llat);   %use same number of points as original data
mlonmin=min(mlondat(:));
mlonmax=max(mlondat(:));
mlon=linspace(mlonmin,mlonmax,llon);
theta=pi/2-mlat*pi/180;    %don't sort since we are making a glon,glat grid out of this...
phi=mlon*pi/180;
[THETA,PHI]=meshgrid(theta,phi);    %because we've arranged the data as lon,lat.
[GLAT,GLON]=geomag2geog(THETA,PHI);
Q=zeros(lt,llon,llat);
E0=zeros(lt,llon,llat);
fprintf('Spatial interpolation step...\n');
for it=1:lt
  Qtmp=squeeze(Qdat(it,:,:));
  Qtmp=interp2(glat,gloncorrected,Qtmp,GLAT(:),GLON(:));    %data are plaid in geographic so do the interpolation in that variable
  Q(it,:,:)=reshape(Qtmp,[1 llon llat]);
  E0tmp=squeeze(E0dat(it,:,:));
  E0tmp=interp2(glat,gloncorrected,E0tmp,GLAT(:),GLON(:));    %data are plaid in geographic so do the interpolation in that variable
  E0(it,:,:)=reshape(E0tmp,[1 llon llat]);  
end


%SET UP TIME VARIABLES
ymd=[2017,03,02];
UTsec=7*3600+time;     %time given in file is the seconds from beginning of hour
UThrs=UTsec/3600;
expdate=cat(2,repmat(ymd,[lt,1]),UThrs,zeros(lt,1),zeros(lt,1));
t=datenum(expdate);


%APPLY SOME SMOOTHING
fprintf('Smoothing longitude...\n');
Qsmooth=zeros(lt,llon,llat);
E0smooth=zeros(lt,llon,llat);
for it=1:lt
    %it
    for ilat=1:llat
        Qtmp=squeeze(Q(it,:,ilat));
        inds=find(isnan(Qtmp));
        Qtmp(inds)=0;
        Qtmp=smooth(Qtmp,4);
        Qsmooth(it,:,ilat)=reshape(Qtmp,[1,llon,1]);
        E0tmp=squeeze(E0(it,:,ilat));
        inds=find(isnan(E0tmp));
        E0tmp(inds)=0;
        E0tmp=smooth(E0tmp,4);
        E0smooth(it,:,ilat)=reshape(E0tmp,[1,llon,1]);
    end
end

fprintf('Smoothing latitude...\n');
for it=1:lt
    %it
    for ilon=1:llon
        Qtmp=smooth(squeeze(Qsmooth(it,ilon,:)),4);
        Qsmooth(it,ilon,:)=reshape(Qtmp,[1,1,llat]);
        E0tmp=smooth(squeeze(E0smooth(it,ilon,:)),4);
        E0smooth(it,ilon,:)=reshape(E0tmp,[1,1,llat]);
    end
end

fprintf('Smoothing time...\n');
for ilat=1:llat
    %ilat
    for ilon=1:llon
        Qtmp=smooth(squeeze(Qsmooth(:,ilon,ilat)),5);
        Qsmooth(:,ilon,ilat)=reshape(Qtmp,[lt,1,1]);
        E0tmp=smooth(squeeze(E0smooth(:,ilon,ilat)),10);    %data are noisier than the total energy flux
        E0smooth(:,ilon,ilat)=reshape(E0tmp,[lt,1,1]);        
    end
end


% %VISUALIZE THE ORIGINAL AND SMOOTHED DATA
% plotdir='./plots/';
% system(['mkdir ',plotdir]);
% figure;
% set(gcf,'PaperPosition',[0 0 8.5 3.5]);
% for it=1:lt
%     clf;
%     subplot(121);
%     imagesc(mlon,mlat,squeeze(Q(it,:,:))');
%     axis xy; axis square;
%     ylabel('mlat.');
%     xlabel('mlon.');
%     title(['Unsmoothed:  ',datestr(datenum(expdate(it,:)))]);
%     colorbar;
%     %cax=caxis;
%     caxis([0 40]);
%     
%     subplot(122);
%     imagesc(mlon,mlat,squeeze(Qsmooth(it,:,:))');
%     axis xy; axis square;
%     ylabel('mlat.');
%     xlabel('mlon.');
%     title(['Smoothed:  ',datestr(datenum(expdate(it,:)))]);
%     colorbar;
%     %caxis(cax);
%     caxis([0 40]);
%     
%     UTsec=expdate(it,4)*3600+expdate(it,5)*60+expdate(it,6);
%     ymd=expdate(it,1:3);
%     filename=datelab(ymd,UTsec);
%     filename=[plotdir,filename,'.png']    
%     
%     print('-dpng',filename,'-r300')
% end
% close all;


%INTERPOLATE/DECIMATE TO 1 SECOND RESOLUTION
%samplerate=1;    %sampling rate in seconds
samplerate=0.5;
startdate=[2017,3,2,7+35/60,0,0];
startt=datenum(startdate);
%tmin=ceil(min(t)*86400)/86400;
tmin=ceil(startt*86400)/86400;
tmax=floor(max(t)*86400)/86400;
outputt=tmin:samplerate/86400:tmax;
outputdate=datevec(outputt);
ltout=numel(outputt);
Qit=zeros(llon,llat,ltout);
E0it=zeros(llon,llat,ltout);
for ilat=1:llat
   for ilon=1:llon
       Qhere=squeeze(Qsmooth(:,ilon,ilat));
       E0here=squeeze(E0smooth(:,ilon,ilat));
       Qi=interp1(t,Qhere,outputt);
       E0i=interp1(t,E0here,outputt);
       Qit(ilon,ilat,:)=reshape(Qi,[1 1 ltout]);
       E0it(ilon,ilat,:)=reshape(E0i,[1 1 ltout]);       
   end
end


% %PLOT THE DECIMATED DATA
% plotdir='./plots_decimated/';
% system(['mkdir ',plotdir]);
% figure;
% for it=1:ltout
%     clf;
%     imagesc(mlon,mlat,squeeze(Qit(:,:,it))');
%     axis xy; axis square;
%     ylabel('mlat.');
%     xlabel('mlon.');
%     title(['Decimated:  ',datestr(datenum(outputdate(it,:)))]);
%     colorbar;
%     cax=caxis;
%     caxis([0 40]);
%     
%     UTsec=outputdate(it,4)*3600+outputdate(it,5)*60+outputdate(it,6);
%     ymd=outputdate(it,1:3);
%     filename=datelab(ymd,UTsec);
%     filename=[plotdir,filename,'.png']
%     
%     print('-dpng',filename,'-r300')
% end
% close all;


%CONVER THE ENERGY TO EV
E0it=max(E0it,1);
E0it=E0it*1e3;


%SAVE THIS DATA TO APPROPRIATE FILES - LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
%FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.  THE EFIELD DATA DO
%NO NEED TO BE SMOOTHED.
filename=[outdir,'simsize.dat'];
fid=fopen(filename,'w');
fwrite(fid,llon,'integer*4');
fwrite(fid,llat,'integer*4');
fclose(fid);
filename=[outdir,'simgrid.dat'];
fid=fopen(filename,'w');
fwrite(fid,mlon,'real*8');
fwrite(fid,mlat,'real*8');
fclose(fid);
for it=1:ltout
    UTsec=outputdate(it,4)*3600+outputdate(it,5)*60+outputdate(it,6);
    ymd=outputdate(it,1:3);
    filename=datelab(ymd,UTsec);
    filename=[outdir,filename,'.dat']
    fid=fopen(filename,'w');
    fwrite(fid,Qit(:,:,it),'real*8');
    fwrite(fid,E0it(:,:,it),'real*8');
    fclose(fid);
end


%ALSO SAVE TO A  MATLAB FILE
save('-v7.3',[outdir,'particles.mat'],'glon','glat','mlon','mlat','Qit','E0it','outputdate');


%RESTORE PATH
rmpath([idlpath,'./restore_idl']);
rmpath ../../script_utils;

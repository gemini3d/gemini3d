addpath ./script_utils;


%OUTPUT FILE LOCATION
outdir='./fields/';
system(['mkdir ',outdir]);


%READ IN FIELD AND POSITION DATA FROM AMISR HDF5 FILE
Exgeog=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Fit2D/Ex_geo'),'double');
Eygeog=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Fit2D/Ey_geo'),'double');
Exgeomagdat=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Fit2D/Ex'),'double');
Eygeomagdat=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Fit2D/Ey'),'double');
Xgeo=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Grid/X_geo'),'double');    % Geo lon
Ygeo=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Grid/Y_geo'),'double');    % Geo lat
[llon,llat,lt]=size(Exgeomagdat);


%TIMING INFORMATION
day=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Time/Day'),'double');
day=mean(day,1);   %average in time
day=day(:);
month=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Time/Month'),'double');
month=mean(month,1);
month=month(:);
year=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Time/Year'),'double');
year=mean(year,1);
year=year(:);
UThrs=cast(hdf5read('20170302.002_lp_1min-cal_2dVEF_001001_NoLonely-geo.h5','Time/dtime'),'double');
UThrs=mean(UThrs,1);
UThrs=UThrs(:);
UTsec=UThrs*3600;
expdate=cat(2,year,month,day,UThrs,zeros(lt,1),zeros(lt,1));    %create a date vector for this dataset
t=datenum(expdate);


%CONVERT POSITIONS TO GEOMAGNETIC USING SAME ALGORITHM AS SIMULATION
glat=Ygeo;
glon=Xgeo;
gloncorrected=360+glon;
% [theta,phi]=geog2geomag(glat,glon);
% mlat=90-theta*180/pi;
% mlon=phi*180/pi;


%SAMPLE THE DATA ON A UNIFORM MLAT,MLON GRID
fprintf('Interpolating...\n')
theta=zeros(llon,llat);
phi=zeros(llon,llat);
for ilat=1:llat
    for ilon=1:llon
        [thetatmp,phitmp]=geog2geomag(glat(ilon,ilat),glon(ilon,ilat));
        theta(ilon,ilat)=thetatmp;
        phi(ilon,ilat)=phitmp;
    end
end
mlatdat=90-theta*180/pi;
mlondat=phi*180/pi;
mlatmin=min(mlatdat(:));
mlatmax=max(mlatdat(:));
mlat=linspace(mlatmin,mlatmax,llat);   %use same number of points as original data
mlonmin=min(mlondat(:));
mlonmax=max(mlondat(:));
mlon=linspace(mlonmin,mlonmax,llon);
[MLAT,MLON]=meshgrid(mlat,mlon);    %mlon goes down the vertical dimension of the matrices in this script
theta=pi/2-mlat*pi/180;    %don't sort since we are making a glon,glat grid out of this...
phi=mlon*pi/180;
[THETA,PHI]=meshgrid(theta,phi);    %because we've arranged the data as lon,lat.
[GLAT,GLON]=geomag2geog(THETA,PHI);
Exgeomag=zeros(llon,llat,lt);
Eygeomag=zeros(llon,llat,lt);
for it=1:lt
  Extmp=squeeze(Exgeomagdat(:,:,it));
  F=TriScatteredInterp(gloncorrected(:),glat(:),Extmp(:));   %the source data are not on a plaid grid here...
  Extmp=F(GLON(:),GLAT(:));
  Exgeomag(:,:,it)=reshape(Extmp,[llon,llat]);
  Eytmp=squeeze(Eygeomagdat(:,:,it));
  F=TriScatteredInterp(gloncorrected(:),glat(:),Eytmp(:));
  Eygeomag(:,:,it)=reshape(F(GLON(:),GLAT(:)),[llon,llat]);
end


% %DEBUG PLOTTING OF ORIGINAL DATA
% plotdir='./plotsE/';
% system(['mkdir ',plotdir]);
% figure;
% set(gcf,'PaperPosition',[0 0 8.5 3.5]);
% for it=1:lt
%     clf;
%     subplot(131);
%     quiver(gloncorrected,glat,Exgeog(:,:,it),Eygeog(:,:,it));
%     xlabel('geographic lon.');
%     ylabel('geographic lat.');
%     title(datestr(datenum(expdate(it,:))));
%     
%     subplot(132);
%     quiver(GLON,GLAT,Exgeomag(:,:,it),Eygeomag(:,:,it));
%     xlabel('geog. lon. (resamp.)');
%     ylabel('geog. lat. (resamp.)');
%     title(datestr(datenum(expdate(it,:))));
% 
%     subplot(133);
% %    quiver(MLON,MLAT,Exgeomag(:,:,it),Eygeomag(:,:,it));
%     quiver(MLON',MLAT',Exgeomag(:,:,it)',Eygeomag(:,:,it)');
%     xlabel('geomagnetic lon.');
%     ylabel('geomagnetic lat.');
%     title(datestr(datenum(expdate(it,:))));    
%     
%     UTsec=expdate(it,4)*3600+expdate(it,5)*60+expdate(it,6);
%     ymd=expdate(it,1:3);
%     filename=datelab(ymd,UTsec);
%     filename=[plotdir,filename,'.png']
%     
%     print('-dpng',filename,'-r300')
% end
% close all;


%NEED TO SAMPLE DATA ON A UNIFORM TEMPORAL GRID FOR THE MODEL
%samplerate=min(diff(t));    %sampling rate in days, this uses the min from the dataset
samplerate=10/86400;    %sampling rate in days
sampleratesec=samplerate*86400;
sampleratesec=round(sampleratesec);
samplerate=sampleratesec/86400;
%outputt=min(t):samplerate:max(t);    %spans data set
%TOIstartdate=[2017,03,02,28200/3600,0,0];    %pick start and end times for the field data times of interest
%TOIenddate=[2017,03,02,28797/3600,0,0];
TOIstartdate=[2017,03,02,7+35/60,0,0];    %pick start and end times for the field data times of interest
TOIenddate=[2017,03,02,8,0,0];
TOIstartt=datenum(TOIstartdate);
TOIendt=datenum(TOIenddate);
outputt=TOIstartt:samplerate:TOIendt+samplerate;
outputdate=datevec(outputt);
ltout=numel(outputt);
Exit=zeros(llon,llat,ltout);
Eyit=zeros(llon,llat,ltout);
for ilat=1:llat
   for ilon=1:llon
       Exhere=squeeze(Exgeomag(ilon,ilat,:));
       inds=find(isnan(Exhere));
       Exhere(inds)=0;    %assumes the simulation is not encapsulating the entire domain of the electric field data
       Eyhere=squeeze(Eygeomag(ilon,ilat,:));
       inds=find(isnan(Eyhere));
       Eyhere(inds)=0;
       Exi=interp1(t,Exhere,outputt);
       Eyi=interp1(t,Eyhere,outputt);
       Exit(ilon,ilat,:)=reshape(Exi,[1 1 ltout]);
       Eyit(ilon,ilat,:)=reshape(Eyi,[1 1 ltout]);       
   end
end


%CONVERT TO V/M
Exit=Exit*1e-3;
Eyit=Eyit*1e-3;


% %PLOT THE INTERPOLATED DATA AS A CHECK
% plotdir='./plotsE_decimated/';
% system(['mkdir ',plotdir]);
% figure;
% for it=1:ltout
%     clf;
%     quiver(MLON',MLAT',Exit(:,:,it)',Eyit(:,:,it)');
%     xlabel('geomagnetic lon.');
%     ylabel('geomagnetic lat.');
%     title(datestr(datenum(outputdate(it,:))));
%     
%     UTsec=outputdate(it,4)*3600+outputdate(it,5)*60+outputdate(it,6);
%     ymd=outputdate(it,1:3);
%     filename=datelab(ymd,UTsec);
%     filename=[plotdir,filename,'.png']
%     
%     print('-dpng',filename,'-r300')
% end


%SAVE THESE DATA TO APPROPRIATE FILES - LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
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

flagdirich=1;   %if 0 data is interpreted as FAC, else we interpret it as potential
for it=1:ltout
    UTsec=outputdate(it,4)*3600+outputdate(it,5)*60+outputdate(it,6);
    ymd=outputdate(it,1:3);
    filename=datelab(ymd,UTsec);
    filename=[outdir,filename,'.dat']
    fid=fopen(filename,'w');

    fwrite(fid,flagdirich,'real*8');
    fwrite(fid,Exit(:,:,it),'real*8');
    fwrite(fid,Eyit(:,:,it),'real*8');
    fwrite(fid,Vminx1it(:,:,it),'real*8');
    fwrite(fid,Vmaxx1it(:,:,it),'real*8');
    fwrite(fid,Vminx2it(:,it),'real*8');
    fwrite(fid,Vmaxx2it(:,it),'real*8');
    fwrite(fid,Vminx3it(:,it),'real*8');
    fwrite(fid,Vmaxx3it(:,it),'real*8');

    fclose(fid);
end


%ALSO CREATE A MATLAB OUTPUT FILE FOR GOOD MEASURE
save([outdir,'fields.mat'],'glon','glat','mlon','mlat','GLAT','GLON','MLAT','MLON','Exit','Eyit','Exgeog','Eygeog','outputdate');


rmpath ./script_utils;

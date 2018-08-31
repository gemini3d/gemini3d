addpath ./script_utils;


%ELECTRIC FIELD DATA
load ./fields/fields.mat;
mlatfields=MLAT;
mlonfields=MLON;
datefields=outputdate;
tfields=datenum(datefields);


%PARTICLE DATA
load ./particles/particles.mat;
mlatparticles=mlat;
mlonparticles=mlon;
dateparticles=outputdate;
tparticles=datenum(dateparticles);  


%INTERPOLATE FIELD DATA ONTO PARTICLE TIME BASIS
llonfields=size(Exit,1);
llatfields=size(Exit,2);
ltpart=numel(tparticles);
Exitpart=zeros(llonfields,llatfields,ltpart);
Eyitpart=zeros(llonfields,llatfields,ltpart);
for ilon=1:llonfields
    for ilat=1:llatfields
      Exitpart(ilon,ilat,:)=interp1(tfields,squeeze(Exit(ilon,ilat,:)),tparticles);
      Eyitpart(ilon,ilat,:)=interp1(tfields,squeeze(Eyit(ilon,ilat,:)),tparticles);      
    end
end


%SIDE BY SIDE PLOTS OF FIELDS AND PARTICLES
plotdir='./plots_synth/';
system(['mkdir ',plotdir]);
figure;
set(gcf,'PaperPosition',[0 0 8.5 3.5]);
% for it=1:ltpart
%     clf;
%     subplot(121);
%     imagesc(mlonparticles,mlatparticles,squeeze(Qit(:,:,it))');
%     axis xy; 
%     axis square;
%     ax=axis;
%     xlabel('mlat.');
%     ylabel('mlon.');    
% 
%     subplot(122);
%     quiver(mlonfields',mlatfields',squeeze(Exitpart(:,:,it))',squeeze(Eyitpart(:,:,it))');
%     xlabel('mlon.');
%     ylabel('mlat.');
%     title(datestr(datenum(dateparticles(it,:))));
%     axis(ax);
%     axis square;
% 
%     UTsec=outputdate(it,4)*3600+outputdate(it,5)*60+outputdate(it,6);
%     ymd=outputdate(it,1:3);
%     filename=datelab(ymd,UTsec);
%     filename=[plotdir,filename,'.png']
%     
%     print('-dpng',filename,'-r300')
% end

% for it=1:ltpart
%     clf;
%     imagesc(mlonparticles,mlatparticles,squeeze(Qit(:,:,it))');
%     colorbar;
%     caxis([0 40])
%     axis xy; 
%     axis square;
%     ax=axis;
%     xlabel('mlat.');
%     ylabel('mlon.');    
%     hold on;   
%     stride=3;
%     quiver(mlonfields(1:stride:end,1:stride:end)',mlatfields(1:stride:end,1:stride:end)',squeeze(Exitpart(1:stride:end,1:stride:end,it))',squeeze(Eyitpart(1:stride:end,1:stride:end,it))','white');
% 
%     UTsec=outputdate(it,4)*3600+outputdate(it,5)*60+outputdate(it,6);
%     ymd=outputdate(it,1:3);
%     filename=datelab(ymd,UTsec);
%     filename=[plotdir,filename,'.png']
%     
%     print('-dpng',filename,'-r300')
% end


for it=1:ltpart
    clf;
    subplot(121);
    imagesc(mlonparticles,mlatparticles,squeeze(Qit(:,:,it))');
    c=colorbar;
    caxis([0 40]);
    ylabel(c,'Q (mW/m^2)');
    axis xy; 
    axis square;
    axis([250 265 65.5 69]);
    ylabel('mlat.');
    xlabel('mlon.');
    title(datestr(datenum(dateparticles(it,:))));
    hold on;   
    %div=divergence(Exitpart(:,:,it)',Eyitpart(:,:,it)');
    %contour(mlonfields',mlatfields',div);
    stride=3;
    quiver(mlonfields(1:stride:end,1:stride:end)',mlatfields(1:stride:end,1:stride:end)',squeeze(Exitpart(1:stride:end,1:stride:end,it))',squeeze(Eyitpart(1:stride:end,1:stride:end,it))','white');

    subplot(122);
    imagesc(mlonparticles,mlatparticles,squeeze(E0it(:,:,it))');
    c=colorbar;
    caxis([0 20e3]);
    ylabel(c,'E_0 (eV)');
    axis xy; 
    axis square;
    axis([250 265 65.5 69]);
    ylabel('mlat.');
    xlabel('mlon.');
    title(datestr(datenum(dateparticles(it,:))));
    hold on;   
    stride=3;
    quiver(mlonfields(1:stride:end,1:stride:end)',mlatfields(1:stride:end,1:stride:end)',squeeze(Exitpart(1:stride:end,1:stride:end,it))',squeeze(Eyitpart(1:stride:end,1:stride:end,it))','white');

    UTsec=outputdate(it,4)*3600+outputdate(it,5)*60+outputdate(it,6);
    ymd=outputdate(it,1:3);
    filename=datelab(ymd,UTsec);
    filename=[plotdir,filename,'.png']
    
    print('-dpng',filename,'-r300')
end    
    

rmpath ./script_utils;

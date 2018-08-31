addpath ../script_utils;

%SIMULATIONS LOCAITONS
simname='tohoku20113D_highres_long/';
basedir='~/zettergmdata/simulations/'
direc=[basedir,simname];
system(['mkdir ',direc,'/magplots']);    %store output plots with the simulation data


%LOAD THE COMPUTED MAGNETIC FIELD DATA
load([direc,'/magfields.mat']);
lt=numel(t);


%INTERPOLATE TO HIGHER SPATIAL RESOLUTION FOR PLOTTING
llonp=500;
llatp=500;
mlonp=linspace(min(mlon),max(mlon),llonp);
mlatp=linspace(min(mlat),max(mlat),llatp);
[MLONP,MLATP]=meshgrid(mlonp,mlatp);
for it=1:lt
    param=interp2(mlon,mlat,squeeze(Brt(:,:,:,it)),MLONP,MLATP);
    Brtp(:,:,:,it)=reshape(param,[1, llonp, llatp]);
    param=interp2(mlon,mlat,squeeze(Bthetat(:,:,:,it)),MLONP,MLATP);
    Bthetatp(:,:,:,it)=reshape(param,[1, llonp, llatp]);
    param=interp2(mlon,mlat,squeeze(Bphit(:,:,:,it)),MLONP,MLATP);
    Bphitp(:,:,:,it)=reshape(param,[1, llonp, llatp]);
end


%SIMULATION META-DATA
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs/config.dat']);


%TABULATE THE SOURCE LOCATION
mlatsrc=mloc(1);
mlonsrc=mloc(2);
thdist=pi/2-mlatsrc*pi/180;    %zenith angle of source location
phidist=mlonsrc*pi/180;


%SETUP FIGURE
figure(1);
set(gcf,'PaperPosition',[0 0 10.5 3]);


%MAKE THE PLOTS AND SAVE TO A FILE
for it=1:lt
    fprintf('Printing magnetic field plots...\n');
    %CREATE A MAP AXIS
    figure(1);
    FS=8;
    
    datehere=simdate_series(it,:);
    ymd=datehere(1:3);
    UTsec=datehere(4)*3600+datehere(5)*60+datehere(6);
    filename=datelab(ymd,UTsec);
    filename=[filename,'.dat'];
    titlestring=datestr(datenum(datehere));
    
    subplot(131);
    cla;
    mlatlimplot=[min(mlat)-0.5,max(mlat)+0.5];
    mlonlimplot=[min(mlon)-0.5,max(mlon)+0.5];
    axesm('MapProjection','Mercator','MapLatLimit',mlatlimplot,'MapLonLimit',mlonlimplot);
    param=squeeze(Brtp(:,:,:,it))*1e9;
    %imagesc(mlon,mlat,param);
    mlatlim=[min(mlatp),max(mlatp)];
    mlonlim=[min(mlonp),max(mlonp)];
    [MLAT,MLON]=meshgrat(mlatlim,mlonlim,size(param));
    pcolorm(MLAT,MLON,param);
    colormap(parula(256));
    set(gca,'FontSize',FS);
    tightmap;
    %axis xy;
    %axis tight;
    caxlim=max(abs(param(:)))
    caxlim=max(caxlim,0.001);
    caxis([-caxlim,caxlim]);
    c=colorbar
    set(c,'FontSize',FS)
    title(['B_r (nT)  ',titlestring])
    xlabel('magnetic long. (deg.)')
    ylabel('magnetic lat. (deg.)')
    hold on;
    ax=axis;
    plotm(mlatsrc,mlonsrc,'r^','MarkerSize',6,'LineWidth',2);
    hold off;
    
    subplot(132);
    cla;
    axesm('MapProjection','Mercator','MapLatLimit',mlatlimplot,'MapLonLimit',mlonlimplot);
    param=squeeze(Bthetatp(:,:,:,it))*1e9;
    %imagesc(mlon,mlat,param);
    pcolorm(MLAT,MLON,param);
    colormap(parula(256));
    set(gca,'FontSize',FS);
    tightmap;
    %axis xy;
    %axis tight;
    caxlim=max(abs(param(:)))
    caxlim=max(caxlim,0.001);
    caxis([-caxlim,caxlim]);
    c=colorbar
    set(c,'FontSize',FS)
    title('B_\theta (nT)')
    xlabel('magnetic long. (deg.)')
    ylabel('magnetic lat. (deg.)')
    hold on;
    ax=axis;
    plotm(mlatsrc,mlonsrc,'r^','MarkerSize',6,'LineWidth',2);
    hold off;
    
    subplot(133);
    cla;
    axesm('MapProjection','Mercator','MapLatLimit',mlatlimplot,'MapLonLimit',mlonlimplot);
    param=squeeze(Bphitp(:,:,:,it))*1e9;
    %imagesc(mlon,mlat,param);
    pcolorm(MLAT,MLON,param);
    colormap(parula(256));
    set(gca,'FontSize',FS);
    tightmap;
    %axis xy;
    %axis tight;
    caxlim=max(abs(param(:)))
    caxlim=max(caxlim,0.001);
    caxis([-caxlim,caxlim]);
    c=colorbar
    set(c,'FontSize',FS)
    title('B_\phi (nT)')
    xlabel('magnetic long. (deg.)')
    ylabel('magnetic lat. (deg.)')
    hold on;
    ax=axis;
    plotm(mlatsrc,mlonsrc,'r^','MarkerSize',6,'LineWidth',2);
    hold off;
    
    subplot(131);
    

    
    %ADD A MAP OF COASTLINES
    if (license('test','Map_Toolbox'))
        load coastlines;
        [thetacoast,phicoast]=geog2geomag(coastlat,coastlon);
        mlatcoast=90-thetacoast*180/pi;
        mloncoast=phicoast*180/pi;
        
        if (360-mlonsrc<20)
            inds=find(mloncoast>180);
            mloncoast(inds)=mloncoast(inds)-360;
        end
        
        subplot(131);
        hold on;
        ax=axis;
        plotm(mlatcoast,mloncoast,'k-','LineWidth',1);
        subplot(132);
        hold on;
        ax=axis;
        plotm(mlatcoast,mloncoast,'k-','LineWidth',1);
        subplot(133);
        hold on;
        ax=axis;
        plotm(mlatcoast,mloncoast,'k-','LineWidth',1);
        hold off;
    end
    axis(ax);
    
    
    print('-dpng',[direc,'/magplots/',filename,'.png'],'-r300');
end

rmpath ../script_utils;

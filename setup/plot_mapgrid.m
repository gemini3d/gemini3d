function ha = plot_mapgrid(xg,flagsource,neuinfo)

narginchk(3,3)
validateattributes(xg, {'struct'}, {'scalar'}, mfilename)
validateattributes(flagsource, {'numeric'}, {'scalar'}, mfilename)
validateattributes(neuinfo, {'struct'}, {'scalar'}, mfilename)
%% INPUT COORDS NEED TO BE CONVERTED TO MAGNETIC
mlon=xg.phi*180/pi;
mlat=90-xg.theta*180/pi;
dmlon=max(mlon(:))-min(mlon(:));


%% ORGANIZE INPUT STRUCTURE
if flagsource ~= 0
    neugridtype=neuinfo.neugridtype;
    sourcelat=neuinfo.sourcelat;
    sourcelong=neuinfo.sourcelong;
    zmin=neuinfo.zmin;
    zmax=neuinfo.zmax;
    if (neugridtype==3)
        xmin=neuinfo.xmin;
        xmax=neuinfo.xmax;
        ymin=neuinfo.ymin;
        ymax=neuinfo.ymax;
    else
        rhomax=neuinfo.rhomax;
    end %if

    [sourcetheta,sourcephi]=geog2geomag(sourcelat,sourcelong);
    sourcemlat=90-sourcetheta*180/pi;
    sourcemlon=sourcephi*180/pi;

    if (360-sourcemlon<dmlon+20)
        sourcemlonplot=sourcemlon-360;
    else
        sourcemlonplot=sourcemlon;
    end
else     %no "epicenter" to track just use mean grid locations
    sourcemlon=mean(mlon(:));
    sourcemlat=mean(mlat(:));
    sourcemlonplot=sourcemlon
end %if


%% SET UP A MAP PLOT IF MAPPING TOOLBOX EXISTS, DETECT WHETHER WE HAVE A 2D OR 3D GRID TO PLOT
flagmapping=license('test','Map_Toolbox');
if (xg.lx(2)==1 || xg.lx(3)==1)    %the grid is 2D
    flag2D=true;
else
    figure;
    hold on;
    flag2D=false;
    if (flagmapping)
        ha=gca;
        axesm('MapProjection','Mercator','MapLatLimit',[-(abs(sourcemlat)+30),abs(sourcemlat)+30],'MapLonLimit',[sourcemlonplot-dmlon/2-10,sourcemlonplot+dmlon/2+10])
        plotfun=@plot3m;
    else
        ha=gca;
        plotfun=@plot3;   %no mapping toolbox so we'll do normal plots - the annoying thing is that these should be lon,lat,alt so all the plots statements have to be switched
    end
end


%% CONVERT INPUT GRID COORDINATES INTO MLAT,MLON,ALT
% Re=6370e3;
% dphi=max(xg.phi(:))-min(xg.phi(:));
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

if (flag2D)    %2D grid, don't use a map to plot and show things in polar coordinates (maybe)
    thref=pi/2;

    h=polar(thref-xg.theta(:,end),xg.r(:,end)/1e3);
    set(h,'LineWidth',LW,'Color',[0 0 0]);
    hold on;
    h=polar(thref-xg.theta(:,1),xg.r(:,1)/1e3);
    set(h,'LineWidth',LW,'Color',[0 0 0]);
    h=polar(thref-xg.theta(1,:),xg.r(1,:)/1e3);
    set(h,'LineWidth',LW,'Color',[0 0 0]);
    h=polar(thref-xg.theta(end,:),xg.r(end,:)/1e3);
    set(h,'LineWidth',LW,'Color',[0 0 0]);
    ha=h;
else    %this is a 3D grid and we can plot it on a map if the user has the appropriate toolbox(es)
    if (flagmapping)    %plots are done lat,lon,alt for plot3m
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
    else    %plotting is done lon,lat,alt with plot3
        h=plotfun(mlon(:,1,1),mlat(:,1,1),alt(:,1,1),'LineWidth',LW);
        plotfun(mlon(:,1,end),mlat(:,1,end),alt(:,1,end),altlinestyle,'LineWidth',LW);
        plotfun(mlon(:,end,1),mlat(:,end,1),alt(:,end,1),'LineWidth',LW);
        h=plotfun(mlon(:,end,end),mlat(:,end,end),alt(:,end,end),altlinestyle,'LineWidth',LW);
        linecolor=h.Color;

        x=squeeze(mlon(1,:,1));
        y=squeeze(mlat(1,:,1));
        z=squeeze(alt(1,:,1));
        plotfun(x,y,z,'LineWidth',LW);

        x=squeeze(mlon(1,:,end));
        y=squeeze(mlat(1,:,end));
        z=squeeze(alt(1,:,end));
        plotfun(x,y,z,altlinestyle,'LineWidth',LW);

        x=squeeze(mlon(end,:,1));
        y=squeeze(mlat(end,:,1));
        z=squeeze(alt(end,:,1));
        plotfun(x,y,z,'LineWidth',LW);

        x=squeeze(mlon(end,:,end));
        y=squeeze(mlat(end,:,end));
        z=squeeze(alt(end,:,end));
        plotfun(x,y,z,altlinestyle,'LineWidth',LW);

        ix3=floor(xg.lx(3)/2);
        x=squeeze(mlon(1,1,1:ix3));
        y=squeeze(mlat(1,1,1:ix3));
        z=squeeze(alt(1,1,1:ix3));
        plotfun(x,y,z,'LineWidth',LW);

        x=squeeze(mlon(1,1,ix3:xg.lx(3)));
        y=squeeze(mlat(1,1,ix3:xg.lx(3)));
        z=squeeze(alt(1,1,ix3:xg.lx(3)));
        plotfun(x,y,z,altlinestyle,'LineWidth',LW);

        x=squeeze(mlon(1,end,1:ix3));
        y=squeeze(mlat(1,end,1:ix3));
        z=squeeze(alt(1,end,1:ix3));
        plotfun(x,y,z,'LineWidth',LW);

        x=squeeze(mlon(1,end,ix3:xg.lx(3)));
        y=squeeze(mlat(1,end,ix3:xg.lx(3)));
        z=squeeze(alt(1,end,ix3:xg.lx(3)));
        plotfun(x,y,z,altlinestyle,'LineWidth',LW);

        x=squeeze(mlon(end,1,1:ix3));
        y=squeeze(mlat(end,1,1:ix3));
        z=squeeze(alt(end,1,1:ix3));
        plotfun(x,y,z,'LineWidth',LW);

        x=squeeze(mlon(end,1,ix3:xg.lx(3)));
        y=squeeze(mlat(end,1,ix3:xg.lx(3)));
        z=squeeze(alt(end,1,ix3:xg.lx(3)));
        plotfun(x,y,z,altlinestyle,'LineWidth',LW);

        x=squeeze(mlon(end,end,1:ix3));
        y=squeeze(mlat(end,end,1:ix3));
        z=squeeze(alt(end,end,1:ix3));
        plotfun(x,y,z,'LineWidth',LW);

        x=squeeze(mlon(end,end,ix3:xg.lx(3)));
        y=squeeze(mlat(end,end,ix3:xg.lx(3)));
        z=squeeze(alt(end,end,ix3:xg.lx(3)));
        plotfun(x,y,z,altlinestyle,'LineWidth',LW);

        xlabel('magnetic latitude (deg.)');
        ylabel('magnetic longitude (deg.)');
        zlabel('altitidue (km)');
    end
    h=gca;     %now go back and make all lines the same color...
    lline=numel(h.Children);
    for iline=1:16    %the last three children are the surface and text label objects
        h.Children(iline).Color=[0 0 0];
    end
end %if


%% PLOT THE NEUTRAL SOURCE LOCATION, IF REQUESTED, right now this assumes axisymmetric
if (flag2D)
    if (flagsource~=0)
        %Create a 2D neutral grid, all distances in km here
        lz=750;
        zn=linspace(zmin,zmax,lz)';
        lrho=750;
        rhomin=0;    %assume origin is included in the neutral simulation
        drho=rhomax-rhomin;
        xn=linspace(-1*drho,drho,lrho);
        rn=zn+6370;
        dtheta=(max(xn(:))-min(xn(:)))/rn(1);
        thetan=linspace(sourcetheta-dtheta/2,sourcetheta+dtheta/2,lrho);
        [THETAn,Rn]=meshgrid(thetan,rn);

        %Now plot the neutral grid
        h=polar(thref-THETAn(:,end),Rn(:,end),altlinestyle);
        set(h,'LineWidth',LW);
        h=polar(thref-THETAn(:,1),Rn(:,1),altlinestyle);
        set(h,'LineWidth',LW);
        h=polar(thref-THETAn(1,:),Rn(1,:),altlinestyle);
        set(h,'LineWidth',LW);
        h=polar(thref-THETAn(end,:),Rn(end,:),altlinestyle);
        set(h,'LineWidth',LW);
    end
else     %full 3D grid
    if (flagsource~=0)    %plot the neutral source and grid
        plotfun(sourcemlat,sourcemlonplot,0,'ro','MarkerSize',16);    %put down a marker where the neutral source will be centered
        lpts=100;                                                     %this is only number of points for the lines to be drawn

        if (neugridtype==2)    %interpret neutral grid as axisymmtric
            %Create a neutral grid, all distance in km here
            rhomin=0;
            zn=linspace(zmin,zmax,lpts);
            rn=6370+zn;    %geocentric distance (in km)

            drho=rhomax-rhomin;                                                  %radius of circle, in kilometers, describing perp. directions of axisymmetric model
            xn=linspace(-1*drho,drho,lpts);                                       %N-S distance spanned by neutral model ("fake" number of grid points used here)
            dthetan=(max(xn(:))-min(xn(:)))/rn(1);                               %equivalent theta coordinates of the neutral mesh (used in the plot of grid)
            thetan=linspace(sourcetheta-dthetan/2,sourcetheta+dthetan/2,lpts);    %theta coordinates of N-S distance specified
            phinhalf1=sourcephi+sqrt((dthetan/2)^2-(thetan-sourcetheta).^2);
            phinhalf2=sourcephi-sqrt((dthetan/2)^2-(thetan-sourcetheta).^2);
            mlatn=90-thetan*180/pi;
            mlonnhalf1=phinhalf1*180/pi;
            mlonnhalf2=phinhalf2*180/pi;
            linemlon=sourcephi*ones(size(zn))*180/pi;

            %Do mlon correction to wrap around
            if (360-sourcemlon<20)
                inds=find(mlonnhalf1>180);
                mlonnhalf1(inds)=mlonnhalf1(inds)-360;
                inds=find(mlonnhalf2>180);
                mlonnhalf2(inds)=mlonnhalf2(inds)-360;
                linemlon=linemlon-360;
            end

            %plot the 3D neutral grid
            hold on;
            zn=zn/altscale;
            plotfun(mlatn,real(mlonnhalf1),zn(1)*ones(size(mlatn)),altlinestyle,'LineWidth',LW);
            plotfun(mlatn,real(mlonnhalf2),zn(1)*ones(size(mlatn)),'LineWidth',LW);
            plotfun(mlatn,real(mlonnhalf1),zn(end)*ones(size(mlatn)),altlinestyle,'LineWidth',LW);
            plotfun(mlatn,real(mlonnhalf2),zn(end)*ones(size(mlatn)),'LineWidth',LW);
            plotfun(min(mlatn)*ones(size(zn)),linemlon,zn,'LineWidth',LW);
            plotfun(max(mlatn)*ones(size(zn)),linemlon,zn,'LineWidth',LW);
            hold off;

            %Make all grid lines the same color
            h=gca;
            lline=numel(h.Children);
            for iline=1:6    %the line objects for each axis are added in a "stack" fashion (last in first out)
                h.Children(iline).Color=linecolor;
            end
        elseif (neugridtype==1)   %we have a 2D Cartesian input neutral grid
            rhomin=0;
            zn=linspace(zmin,zmax,lpts);
            drho=rhomax-rhomin;
            xn=linspace(-drho,drho,lpts);
            yn=xn;              %arbitrarily let the y-extent be the same as x (our defs. of x,y left-handed here)
            rn=zn+6370;         %convert altitude to geocentric distance

            dtheta=(max(xn(:))-min(xn(:)))/rn(1);    %equivalent theta coordinates of the neutral mesh
            dphi=(max(yn(:))-min(yn(:)))/rn(1)/sin(sourcetheta);
            thetan=linspace(sourcetheta-dtheta/2,sourcetheta+dtheta/2,lpts);
            phin=linspace(sourcephi-dphi/2,sourcephi+dphi/2,lpts);
            [THETAn,PHIn,Rn]=meshgrid(thetan,phin,rn);

            MLATn=90-THETAn*180/pi;
            MLONn=PHIn*180/pi;
            Zn=(Rn-6370);
            Zn=Zn/altscale;     %must scale this if it is going on a map plot

            hold on;
            plotfun(MLATn(:,end,1),MLONn(:,end,1),Zn(:,end,1),'LineWidth',LW);
            h=plotfun(MLATn(:,end,end),MLONn(:,end,end),Zn(:,end,end),'LineWidth',LW);
            linecolor=h.Color;    %grab the color of the second line
            plotfun(MLATn(:,1,1),MLONn(:,1,1),Zn(:,1,1),'LineWidth',LW);
            plotfun(MLATn(:,1,end),MLONn(:,1,end),Zn(:,1,end),'LineWidth',LW);
            plotfun(squeeze(MLATn(1,:,1)),squeeze(MLONn(1,:,1)),squeeze(Zn(1,:,1)),'LineWidth',LW);
            plotfun(squeeze(MLATn(1,:,end)),squeeze(MLONn(1,:,end)),squeeze(Zn(1,:,end)),'LineWidth',LW);
            plotfun(squeeze(MLATn(end,:,1)),squeeze(MLONn(end,:,1)),squeeze(Zn(end,:,1)),altlinestyle,'LineWidth',LW);
            plotfun(squeeze(MLATn(end,:,end)),squeeze(MLONn(end,:,end)),squeeze(Zn(end,:,end)),altlinestyle,'LineWidth',LW);
            plotfun(squeeze(MLATn(1,1,:)),squeeze(MLONn(1,1,:)),squeeze(Zn(1,1,:)),'LineWidth',LW);
            plotfun(squeeze(MLATn(1,end,:)),squeeze(MLONn(1,end,:)),squeeze(Zn(1,end,:)),'LineWidth',LW);
            plotfun(squeeze(MLATn(end,1,:)),squeeze(MLONn(end,1,:)),squeeze(Zn(end,1,:)),altlinestyle,'LineWidth',LW);
            plotfun(squeeze(MLATn(end,end,:)),squeeze(MLONn(end,end,:)),squeeze(Zn(end,end,:)),altlinestyle,'LineWidth',LW);

            %Make all grid lines the same color
            h=gca;
            lline=numel(h.Children);
            for iline=1:12    %the line objects for each axis are added in a "stack" fashion (last in first out)
                h.Children(iline).Color=linecolor;
            end
        elseif (neugridtype==3)  %3D Cartesian input neutral grid
            zn=linspace(zmin,zmax,lpts);
            xn=linspace(xmin,xmax,lpts);
            yn=linspace(ymin,ymax,lpts);
            rn=zn+6370;         %convert altitude to geocentric distance

            dtheta=(max(xn(:))-min(xn(:)))/rn(1);    %equivalent theta coordinates of the neutral mesh
            dphi=(max(yn(:))-min(yn(:)))/rn(1)/sin(sourcetheta);
            thetan=linspace(sourcetheta-dtheta/2,sourcetheta+dtheta/2,lpts);
            phin=linspace(sourcephi-dphi/2,sourcephi+dphi/2,lpts);
            [THETAn,PHIn,Rn]=meshgrid(thetan,phin,rn);

            MLATn=90-THETAn*180/pi;
            MLONn=PHIn*180/pi;
            Zn=(Rn-6370);
            Zn=Zn/altscale;     %must scale this if it is going on a map plot

            hold on;
            plotfun(MLATn(:,end,1),MLONn(:,end,1),Zn(:,end,1),'LineWidth',LW);
            h=plotfun(MLATn(:,end,end),MLONn(:,end,end),Zn(:,end,end),'LineWidth',LW);
            linecolor=h.Color;    %grab the color of the second line
            plotfun(MLATn(:,1,1),MLONn(:,1,1),Zn(:,1,1),'LineWidth',LW);
            plotfun(MLATn(:,1,end),MLONn(:,1,end),Zn(:,1,end),'LineWidth',LW);
            plotfun(squeeze(MLATn(1,:,1)),squeeze(MLONn(1,:,1)),squeeze(Zn(1,:,1)),'LineWidth',LW);
            plotfun(squeeze(MLATn(1,:,end)),squeeze(MLONn(1,:,end)),squeeze(Zn(1,:,end)),'LineWidth',LW);
            plotfun(squeeze(MLATn(end,:,1)),squeeze(MLONn(end,:,1)),squeeze(Zn(end,:,1)),altlinestyle,'LineWidth',LW);
            plotfun(squeeze(MLATn(end,:,end)),squeeze(MLONn(end,:,end)),squeeze(Zn(end,:,end)),altlinestyle,'LineWidth',LW);
            plotfun(squeeze(MLATn(1,1,:)),squeeze(MLONn(1,1,:)),squeeze(Zn(1,1,:)),'LineWidth',LW);
            plotfun(squeeze(MLATn(1,end,:)),squeeze(MLONn(1,end,:)),squeeze(Zn(1,end,:)),'LineWidth',LW);
            plotfun(squeeze(MLATn(end,1,:)),squeeze(MLONn(end,1,:)),squeeze(Zn(end,1,:)),altlinestyle,'LineWidth',LW);
            plotfun(squeeze(MLATn(end,end,:)),squeeze(MLONn(end,end,:)),squeeze(Zn(end,end,:)),altlinestyle,'LineWidth',LW);

            %Make all grid lines the same color
            h=gca;
            lline=numel(h.Children);
            for iline=1:12    %the line objects for each axis are added in a "stack" fashion (last in first out)
                h.Children(iline).Color=linecolor;
            end
        end %if
    end %if
end %if


%% CLEAN UP FONT SIZES, ETC.
FS=12;
set(gca,'FontSize',FS);
grid on;
set(gca,'LineWidth',LW-0.5)


%% ADD MAPPED COASTLINES TO PLOT - DONE LAST TO MAKE SOME OF THE GRID LINE PLOTTING CODE CLEANER
if (flag2D)    %for 2D just plot Earth's surface for reference
    terr_theta=linspace(0,2*pi,100);
    terr_r=6370*ones(size(terr_theta));
    h=polar(terr_theta,terr_r,'c:');
    set(h,'LineWidth',LW);
else
    if (flagmapping)
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
end %if



%% OLD/EXTRA CODE THAT I'M LOATH TO PART WITH (MZ)
% %MAKE A MOVIE OF THE GRID ROTATING
% direc='~/Downloads/gridplot/';
% makedir(direc)
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

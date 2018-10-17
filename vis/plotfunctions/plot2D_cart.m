function h=plot2D_cart(ymd,UTsec,xg,parm,parmlbl,caxlims,sourceloc,ha)

narginchk(4,8)
validateattributes(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 1)
validateattributes(UTsec, {'numeric'}, {'scalar'}, mfilename, 'UTC second', 2)
validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 3)
validateattributes(parm, {'numeric'}, {'2d'}, mfilename)
if nargin<5, parmlbl=''; end
validateattributes(parmlbl, {'char'}, {'vector'}, mfilename, 'parameter label', 5)
if nargin<6
  caxlims=[];
else
  validateattributes(caxlims, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'plot intensity (min, max)', 6)
end
if nargin<7  || isempty(sourceloc) % leave || for validate
  sourceloc = [];
else
  validateattributes(sourceloc, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'source magnetic coordinates', 7)
end
if nargin<8 || isempty(ha)
  clf, h=gcf; ha=gca; 
else
  h = ancestor(ha,'figure');
end

%set(h,'PaperPosition',[0 0 11 4.5]);


%REORGANIZE INPUT
dmy = flip(ymd);
t=UTsec/3600;


%SOURCE LOCATION
if ~isempty(sourceloc)
  sourcemlat=sourceloc(1);
  sourcemlon=sourceloc(2);
else
  sourcemlat=[];
  sourcemlon=[];
end


%SIZE OF SIMULATION
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
inds1=3:lx1+2;
inds2=3:lx2+2;
inds3=3:lx3+2;
Re=6370e3;


%JUST PICK AN X3 LOCATION FOR THE MERIDIONAL SLICE PLOT, AND AN ALTITUDE FOR THE LAT./LON. SLICE
ix3=floor(lx3/2);
altref=300;


%SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
meantheta=mean(xg.theta(:));
meanphi=mean(xg.phi(:));
y=-1*(xg.theta-meantheta);   %this is a mag colat. coordinate and is only used for defining grid in linspaces below, runs backward from north distance, hence the negative sign
%x=(xg.phi-meanphi);       %mag. lon coordinate, pos. eastward
x=xg.x2/Re/sin(meantheta);
z=xg.alt/1e3;
lxp=500;
lyp=500;
lzp=500;
minx=min(x(:));
maxx=max(x(:));
miny=min(y(:));
maxy=max(y(:));
minz=min(z(:));
maxz=max(z(:));
xp=linspace(minx,maxx,lxp);     %eastward distance (rads.)
yp=linspace(miny,maxy,lyp);     %should be interpreted as northward distance (in rads.).  Irrespective of ordering of xg.theta, this will be monotonic increasing!!!
zp=linspace(minz,maxz,lzp)';     %altitude (meters)


%INTERPOLATE ONTO PLOTTING GRID
if (xg.lx(3)==1)     %alt./lon. slice
  [X,Z]=meshgrid(xp,zp*1e3);    %meridional meshgrid, this defines the grid for plotting
  x1plot=Z(:);   %upward distance
  x2plot=X(:)*Re*sin(meantheta);     %eastward distance

  parmtmp=parm(:,:);
  parmp=interp2(xg.x2(inds2),xg.x1(inds1),parmtmp,x2plot,x1plot);
  parmp=reshape(parmp,[lzp,lxp]);    %slice expects the first dim. to be "y" ("z" in the 2D case)
elseif (xg.lx(2)==1)     %alt./lat. slice
  [Y3,Z3]=meshgrid(yp,zp*1e3);

  x1plot=Z3(:);   %upward distance
  x3plot=Y3(:)*Re;     %northward distance;

  ix2=floor(lx2/2);
  parmtmp=parm(:,:);     %so north dist, east dist., alt.
  parmp3=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
  parmp3=reshape(parmp3,[lzp,lyp]);    %slice expects the first dim. to be "y"
end


%CONVERT ANGULAR COORDINATES TO MLAT,MLON
if (xg.lx(2)==1)
  yp=yp*Re;
  [yp,inds]=sort(yp);
  %parmp2=parmp2(inds,:,:);
  parmp3=parmp3(:,inds);
end 

if (xg.lx(3)==1)
  %xp=(xp+meanphi)*180/pi;
  xp=xp*Re*sin(meantheta);    %eastward ground distance
  [xp,inds]=sort(xp);
  parmp=parmp(:,inds,:);
  %parmp2=parmp2(:,inds,:);
end

%COMPUTE SOME BOUNDS FOR THE PLOTTING
minxp=min(xp(:));
maxxp=max(xp(:));
minyp=min(yp(:));
maxyp=max(yp(:));
minzp=min(zp(:));
maxzp=max(zp(:));


%GLOBAL PLOT COMMANDS
FS=8;

%MAKE THE PLOT!
if (xg.lx(3)==1)
  hi=imagesc(ha,xp/1e3,zp,parmp);
  hold(ha, 'on')
  plot(ha,[minxp,maxxp],[altref,altref],'w--','LineWidth',2);
  if ~isempty(sourcemlat)
    plot(ha,sourcemlat,0,'r^','MarkerSize',12,'LineWidth',2);
  end
  hold(ha, 'off')
  try % not yet in Octave 4.4.1
    set(hi,'alphadata',~isnan(parmp))
  end
  set(ha,'FontSize',FS)
  axis(ha, 'xy')
  axis(ha, 'square')
  colormap(ha, parula(256))
  if ~isempty(caxlims)
    caxis(ha, caxlims)
  end
  c=colorbar(ha);
  xlabel(c,parmlbl)
  xlabel(ha, 'eastward dist. (km)')
  ylabel(ha, 'altitude (km)')
elseif (xg.lx(2)==1)
  hi=imagesc(ha,yp/1e3,zp,parmp3);
  hold(ha, 'on')
  %plot([minyp,maxyp],[altref,altref],'w--','LineWidth',2);
  if (~isempty(sourcemlat))
    plot(ha, sourcemlat,0,'r^','MarkerSize',12,'LineWidth',2);
  end
  hold(ha, 'off')
  set(hi,'alphadata',~isnan(parmp3));
  set(ha,'FontSize',FS);
  axis(ha, 'xy')
  axis(ha, 'square')
  colormap(ha, parula(256));
  if ~isempty(caxlims)
    caxis(ha, caxlims)
  end
  c=colorbar(ha);
  xlabel(c,parmlbl)
  xlabel(ha, 'northward dist. (km)')
  ylabel(ha, 'altitude (km)')
end


%CONSTRUCT A STRING FOR THE TIME AND DATE
UThrs=floor(t);
UTmin=floor((t-UThrs)*60);
UTsec=floor((t-UThrs-UTmin/60)*3600);
UThrsstr=num2str(UThrs);
UTminstr=num2str(UTmin);
if (numel(UTminstr)==1)
  UTminstr=['0',UTminstr];
end
UTsecstr=num2str(UTsec);
if (numel(UTsecstr)==1)
  UTsecstr=['0',UTsecstr];
end

timestr=[UThrsstr,':',UTminstr,':',UTsecstr];
%strval=sprintf('%s \n %s',[num2str(dmy(1)),'/',num2str(dmy(2)),'/',num2str(dmy(3))], ...
%    [num2str(t),' UT']);
strval=sprintf('%s \n %s',[num2str(dmy(2)),'/',num2str(dmy(1)),'/',num2str(dmy(3))], ...
    [timestr,' UT']);
%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',18,'Color',[0.66 0.66 0.66],'FontWeight','bold');
%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',16,'Color',[0.5 0.5 0.5],'FontWeight','bold');
title(strval)

end
